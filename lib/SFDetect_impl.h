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

#ifndef INCLUDED_SFDETECTOR_SFDETECT_IMPL_H
#define INCLUDED_SFDETECTOR_SFDETECT_IMPL_H

#include <SFDetector/SFDetect.h>
#include <liquid/liquid.h>
#include <string>
#include <vector>
#include <fstream>
#include <volk/volk.h>


namespace gr {
  namespace SFDetector {

    class SFDetect_impl : public SFDetect
    {
     private:
      // Nothing to declare in this block.
      uint32_t d_samp_per_second;
      uint32_t d_freq_off;
      uint32_t d_bandwidth;
      uint32_t d_samp_per_symbol;
      double d_symb_per_second;
      uint32_t d_total_bins;
      double d_dt;
      
      
      std::vector<gr_complex> d_downchirp;        ///< The complex ideal downchirp.
      std::vector<float>      d_downchirp_ifreq;  ///< The instantaneous frequency of the ideal downchirp.

      std::vector<gr_complex> d_upchirp;          ///< The complex ideal upchirp.
      std::vector<float>      d_upchirp_ifreq;    ///< The instantaneous frequency of the ideal upchirp.
      std::vector<float>      d_upchirp_ifreq_v;  ///< The instantaneous frequency of the ideal upchirp.
      std::vector<gr_complex> d_fft;              ///< Vector containing the FFT resuls.
      std::vector<gr_complex> d_mult_hf;          ///< Vector containing the FFT decimation.
      std::vector<gr_complex> d_tmp;              ///< Vector containing the FFT decimation.

      fftplan d_q;                                ///< The LiquidDSP::FFT_Plan.
      fftplan d_qr;                               ///< The LiquidDSP::FFT_Plan in reverse.
 
      void build_ideal_chirps(void);
       /**
                 *  \brief  Calculate the instantaneous frequency for the given complex symbol.
                 *
                 *  \param  in_samples
                 *          The complex array to calculate the instantaneous frequency for.
                 *  \param  out_ifreq
                 *          The output `float` array containing the instantaneous frequency.
                 *  \param  window
                 *          The size of said arrays.
                 */
     inline void instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window);

     uint32_t get_shift_fft(const gr_complex *samples);

     public:
      SFDetect_impl(float samp_rate, float center_freq, std::vector<float> channel_list, uint32_t bandwidth);
      ~SFDetect_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace SFDetector
} // namespace gr

#endif /* INCLUDED_SFDETECTOR_SFDETECT_IMPL_H */

