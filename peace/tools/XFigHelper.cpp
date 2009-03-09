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

XFigHelper::XFigHelper() {
    // Nothing to be done for now.
}

XFigHelper::~XFigHelper() {
    // Nothing else to be done.
}

bool
XFigHelper::setOutput(const std::string &fileName,
                      const bool genCustomColors) {
    os.open(fileName.c_str());
    if (!os.good()) {
        return false;
    }
    // Dump the XFIG header.
    dumpHeader(genCustomColors);
    // Everything went well.
    return true;
}


void
XFigHelper::dumpHeader(const bool genCustomColors) {
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

    if (genCustomColors) {
        generateColorTable();
    }
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
    const int FontHeight = 1200 * fontSize / 72;
    const int BaseLine   = (int) (FontHeight * 0.75);
    
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

void
XFigHelper::generateColorTable() {
    // Dump custom color codes for further use.
    int colorCode = 32;
    for(int red = 0; (red < 256); red += 85) {
        for(int blue = 0; (blue < 256); blue += 85) {
            for(int green = 0; (green < 256); green += 85) {
                char colorString[16];
                sprintf(colorString, "#%02x%02x%02x", red, green, blue);
                os << "0 " << colorCode << " "
                   << colorString << "\n";
                colorCode++;
            }
        }
    }    
}

#endif
