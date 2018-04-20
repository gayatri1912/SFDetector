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


#ifndef INCLUDED_SFDETECTOR_SFDETECT_H
#define INCLUDED_SFDETECTOR_SFDETECT_H

#include <SFDetector/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace SFDetector {

    /*!
     * \brief <+description of block+>
     * \ingroup SFDetector
     *
     */
    class SFDETECTOR_API SFDetect : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<SFDetect> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of SFDetector::SFDetect.
       *
       * To avoid accidental use of raw pointers, SFDetector::SFDetect's
       * constructor is in a private implementation
       * class. SFDetector::SFDetect::make is the public interface for
       * creating new instances.
       */
      static sptr make(float samp_rate, float center_freq, std::vector<float> channel_list, uint32_t bandwidth);
    };

  } // namespace SFDetector
} // namespace gr

#endif /* INCLUDED_SFDETECTOR_SFDETECT_H */

