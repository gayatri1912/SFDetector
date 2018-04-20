/* -*- c++ -*- */

#define SFDETECTOR_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "SFDetector_swig_doc.i"

%{
#include "SFDetector/SFDetect.h"
%}


%include "SFDetector/SFDetect.h"
GR_SWIG_BLOCK_MAGIC2(SFDetector, SFDetect);
