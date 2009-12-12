#ifndef COMMON_H
#define COMMON_H

//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

/** \file Common.h

    \brief Common definitions and constants used by various tools.

    This file contains a collection of common constants and macros
    that are shared and used by multiple tools.
*/

// Font code for COURIER font
#define COURIER 14

// The size of font (in points) to use
#define FONT_SIZE 6

// Some color codes to make code a bit more readable
#define BLACK  0
#define BLUE   1
#define GREEN  2
#define CYAN   3
#define RED    4
#define PINK   5
#define YELLOW 6
#define WHITE  7

// Color code offset from where user colors start
#define USER_COLOR_START 32

#define CHECK_ARGS(tool, condition, errorMsg)	\
  if (condition) {				\
    std::cout << errorMsg;			\
    Tool::showUsage(tool, ap);			\
    return 1;					\
  }

#endif
