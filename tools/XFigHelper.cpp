#ifndef XFIG_HELPER_CPP
#define XFIG_HELPER_CPP

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "XFigHelper.h"

XFigHelper::XFigHelper(std::ostream &outputStream) : os(outputStream) {
    // Dump the XFIG header.
    dumpHeader();
}

XFigHelper::~XFigHelper() {
    // Nothing else to be done.
}

void
XFigHelper::dumpHeader() const {
    os << "#FIG 3.2  Produced by PEACE Tools 0.1\n"
       << "Landscape\n"
       << "Center\n"
       << "Inches\n"
       << "Letter\n"
       << "100.00\n"
       << "Single\n"
       << "-2\n"
       << "1200 2\n"
       << "0 32 #cccccc\n";
}

void
XFigHelper::drawLine(const int x1, const int y1,
                     const int x2, const int y2,
                     const int colorCode, const int level) {
    os  << "2 1 0 1 " << colorCode
        << " 7 " << level << " -1 -1 0.000 0 0 -1 1 0 2\n"
        << "\t 1 1 0.00 0.00 0.00\n"
        << "\t  " << x1 << " " << y1
        << " "    << x2  << " " << y2
        << std::endl; 
}

int
XFigHelper::drawText(const std::string& text, const int x, const int y,
                     const int fontCode, const int fontSize,
                     const int colorCode, const int level) {
    const double FontHeight = (1200.0 / 72.0 * fontSize);
    const int BaseLine      = (int) (FontHeight * 0.75);
    
    os << "4 0 " << colorCode << " " << level << " -1 " << fontCode
       << " "    << fontSize  << " 0.0000 4 0 0 "
       << x << " " << y + BaseLine << " " << text << "\\001\n";
    // Return font height.
    return FontHeight;
}

void
XFigHelper::drawRect(int x, int y, int width, int height, int colorCode) {
    os << "2 2 0 1 " << colorCode << " " << colorCode
       << " 50 -1 20 0.000 0 0 -1 0 0 5\n"
       << "\t" << x           << " " << y
       << " "  << (x + width) << " " << y
       << " "  << (x + width) << " " << (y + height)
       << " "  << x           << " " << (y + height)
       << "\t" << x           << " " << y
       << std::endl;
}

#endif
