/* -*- c++ -*- */
/* 
 * Copyright 2018 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "SFDetect_impl.h"
#include <liquid/liquid.h>
#include <numeric>
#include <algorithm>
#include <gnuradio/expj.h>

namespace gr {
  namespace SFDetector {

    SFDetect::sptr
    SFDetect::make(float samp_rate, float center_freq, std::vector<float> channel_list, uint32_t bandwidth)
    {
      return gnuradio::get_initial_sptr
        (new SFDetect_impl(samp_rate, center_freq, channel_list, bandwidth));
    }

    /*
     * The private constructor
     */
    SFDetect_impl::SFDetect_impl(float samp_rate, float center_freq, std::vector<float> channel_list, uint32_t bandwidth)
      : gr::block("SFDetect",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(int)))
    {
       d_samp_per_second = samp_rate;
       d_freq_off = channel_list[0]-center_freq;
       d_bandwidth = bandwidth;
       d_dt = 1.0f/ d_samp_per_second;
    }

    /*
SFDetect_implSFDetect_implSFDetect_impl     * Our virtual destructor.
     */
    SFDetect_impl::~SFDetect_impl()
    {
    }
    
    inline void SFDetect_impl::instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window) {
      if (window < 2u) {
            std::cerr << "[LoRa Decoder] WARNING : window size < 2 !" << std::endl;
             return;
           }

            /* instantaneous_phase */
       for (uint32_t i = 1u; i < window; i++) {
 	       const float iphase_1 = std::arg(in_samples[i - 1]);
                   float iphase_2 = std::arg(in_samples[i]);

                // Unwrapped loops from liquid_unwrap_phase
                while ( (iphase_2 - iphase_1) >  M_PI ) iphase_2 -= 2.0f*M_PI;
                while ( (iphase_2 - iphase_1) < -M_PI ) iphase_2 += 2.0f*M_PI;

                out_ifreq[i - 1] = iphase_2 - iphase_1;
            }

            // Make sure there is no strong gradient if this value is accessed by mistake
            out_ifreq[window - 1] = out_ifreq[window - 2];
        }


    void SFDetect_impl::build_ideal_chirps(void) {
            d_downchirp.resize(d_samp_per_symbol);
            d_upchirp.resize(d_samp_per_symbol);
            d_downchirp_ifreq.resize(d_samp_per_symbol);
            d_upchirp_ifreq.resize(d_samp_per_symbol);
            d_upchirp_ifreq_v.resize(d_samp_per_symbol*3);
            gr_complex tmp[d_samp_per_symbol*3];

            const double T       = -0.5 * d_bandwidth * d_symb_per_second;
            const double f0      = (d_bandwidth / 2.0);
            const double pre_dir = 2.0 * M_PI;
            double t;
            gr_complex cmx       = gr_complex(1.0f, 1.0f);

            for (uint32_t i = 0u; i < d_samp_per_symbol; i++) {
                // Width in number of samples = samples_per_symbol
                // See https://en.wikipedia.org/wiki/Chirp#Linear
                t = d_dt * i;
                d_downchirp[i] = cmx * gr_expj(pre_dir * t * (f0 + T * t));
                d_upchirp[i]   = cmx * gr_expj(pre_dir * t * (f0 + T * t) * -1.0f);
            }

            // Store instantaneous frequency
            instantaneous_frequency(&d_downchirp[0], &d_downchirp_ifreq[0], d_samp_per_symbol);
            instantaneous_frequency(&d_upchirp[0],   &d_upchirp_ifreq[0],   d_samp_per_symbol);

           // samples_to_file("/tmp/downchirp", &d_downchirp[0], d_downchirp.size(), sizeof(gr_complex));
            //samples_to_file("/tmp/upchirp",   &d_upchirp[0],   d_upchirp.size(),   sizeof(gr_complex));

            // Upchirp sequence
            memcpy(tmp, &d_upchirp[0], sizeof(gr_complex) * d_samp_per_symbol);
            memcpy(tmp+d_samp_per_symbol, &d_upchirp[0], sizeof(gr_complex) * d_samp_per_symbol);
            memcpy(tmp+d_samp_per_symbol*2, &d_upchirp[0], sizeof(gr_complex) * d_samp_per_symbol);
            instantaneous_frequency(tmp, &d_upchirp_ifreq_v[0], d_samp_per_symbol*3);
        }

        uint32_t SFDetect_impl::get_shift_fft(const gr_complex *samples) {
            float fft_mag[d_total_bins];

            //samples_to_file("/tmp/data", &samples[0], d_samp_per_symbol, sizeof(gr_complex));

            // Multiply with ideal downchirp
            for (uint32_t i = 0u; i < d_samp_per_symbol; i++) {
                d_mult_hf[i] = samples[i] * d_downchirp[i];
            }

            //samples_to_file("/tmp/mult", &d_mult_hf[0], d_samp_per_symbol, sizeof(gr_complex));

            // Perform FFT
            fft_execute(d_q);

            // Decimate. Note: assumes fft size is multiple of decimation factor and number of bins is even
            // This decimation should be identical to numpy's approach
            const uint32_t N = d_total_bins;
            memcpy(&d_tmp[0],               &d_fft[0],                                     (N + 1u) / 2u * sizeof(gr_complex));
            memcpy(&d_tmp[ (N + 1u) / 2u ], &d_fft[d_samp_per_symbol - (N / 2u)],        N / 2u * sizeof(gr_complex));
            d_tmp[N / 2u] += d_fft[N / 2u];

            // Get magnitude
            for (uint32_t i = 0u; i < d_total_bins; i++) {
                fft_mag[i] = std::abs(d_tmp[i]);
            }

            //samples_to_file("/tmp/fft", &d_tmp[0], d_total_bins, sizeof(gr_complex));

            fft_execute(d_qr); // For debugging
            //samples_to_file("/tmp/resampled", &d_mult_hf[0], d_total_bins, sizeof(gr_complex));

            // Return argmax here
            return (std::max_element(fft_mag, fft_mag + d_total_bins) - fft_mag);
        }
/* unsigned short
    demod_impl::argmax(gr_complex *fft_result,
                       bool update_squelch)
    {
      float magsq   = pow(real(fft_result[0]), 2) + pow(imag(fft_result[0]), 2);
      float max_val = magsq;
      unsigned short   max_idx = 0;


      for (unsigned short i = 0; i < d_fft_size; i++)
      {
        magsq = pow(real(fft_result[i]), 2) + pow(imag(fft_result[i]), 2);
        if (magsq > max_val)
        {
          max_idx = i;
          max_val = magsq;
        }
      }

      if (update_squelch)
      {
        d_power = max_val;
        d_squelched = (d_power > d_threshold) ? false : true;
      }

      return max_idx;
    }
*/
    void
    SFDetect_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
       ninput_items_required[0] = noutput_items * d_samp_per_second;
    }

    int
    SFDetect_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      int *out = (int *) output_items[0];
      
      for(uint32_t sf=7u ; sf <= 12u ; sf++)
      {
	      int k = 0;
	      uint32_t bin_idx[6];
	      d_symb_per_second = (double) d_bandwidth / (1u << sf);
	      d_samp_per_symbol = (uint32_t) (d_samp_per_second/d_symb_per_second);
	      d_total_bins = (uint32_t)(1u << sf);
	      build_ideal_chirps();
	      bin_idx[k] = get_shift_fft(in);
	      std::cout << "SF " << sf << "Bin index " << bin_idx[k] << std::endl;
	      k++;
      }
      
 
      // Do <+signal processing+>
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace SFDetector */
} /* namespace gr */

