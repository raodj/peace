#ifndef ARG_PARSER_CPP
#define ARG_PARSER_CPP

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

#include "ArgParser.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <sstream>

#define MAXIMUM_ARGUMENTS 99

ArgParser::ArgParser() {
}

ArgParser::ArgParser(const ArgRecord validArguments[]) {
    addValidArguments(validArguments);
}

void 
ArgParser::addValidArguments(const ArgRecord validArguments[]) {
    int scannedArgumentCount = 0;
    // Add the supplied valid arguments to our vector of valid
    // arguments.
    while (true) {
        if (validArguments[scannedArgumentCount].type == INVALID) {
            break;
        } else {
            argRecords.push_back(validArguments[scannedArgumentCount]);
            scannedArgumentCount++;
            // Do a simple sanity check
            if (argRecords.size() > MAXIMUM_ARGUMENTS) {
                std::cerr << "Error: There seem to be more than "
                          << MAXIMUM_ARGUMENTS
                          << " arguments specified" << std::endl
                          << "Please check your code for the "
                          << "{\"\", \"\", NULL, ArgType::INVALID}"
                          << " marker or increase MAXIMUM_ARGUMENTS"
                          << std::endl;
                exit(-1);
            }
        }
    }
}

void 
ArgParser::parseArguments(int& argc, char* argv[], bool caxoe) {
    // This loop cycles through the arguments.
    for (int argument = 1; argument < argc; argument++) {
        // This loop compares the arguments passed in with those we're
        // checking them against
        const int NumArgs = argRecords.size();
        bool foundMatch   = false;
        for (int match = 0; ((match < NumArgs) && (argument < argc)); match++){
            ArgRecord& argRec = argRecords[match];
            if (argRec.command != argv[argument]) {
                // This isn't the correct/valid argument. We should
                // continue trying until we find a valid one.
                continue;
            }
            foundMatch = true;
            // For non-boolean options (which don't need a value)
            // ensure we have a value specified for the argument.
            if ((argRec.type != BOOLEAN) && (argRec.type != STRING_LIST)) {
                if ((argument + 1 == argc) ||
                    ((argv[argument + 1][0] == '-') &&
                     (argv[argument + 1][1] == '-'))) {
                    std::cerr << "Parameter for " << argv[argument]
                              << " is missing" << std::endl;
                    exit(-1);
                }
            }
            
            switch (argRec.type) {
            case BOOLEAN:
                removeArgument(argument, argc, argv);
                *(reinterpret_cast<bool*>(argRec.data)) = true;
                break;
                
            case INTEGER:
                removeArgument(argument, argc, argv);
                *(reinterpret_cast<int*>(argRec.data)) = atoi(argv[argument]);
                removeArgument(argument, argc, argv);
                break;

            case UNSIGNED_INT:
                removeArgument(argument, argc, argv);
                sscanf(argv[argument], "%u",
                       reinterpret_cast<unsigned int*>(argRec.data));
                removeArgument(argument, argc, argv);
                break;
               
            case STRING:
                removeArgument(argument, argc, argv);
                *(reinterpret_cast<std::string*>(argRec.data)) = argv[argument];
                removeArgument(argument, argc, argv);
                break;

            case FLOAT:
                removeArgument(argument, argc, argv);
                *(reinterpret_cast<float*>(argRec.data))= atof(argv[argument]);
                removeArgument(argument, argc, argv);
                break;
                
            case DOUBLE:
                removeArgument(argument, argc, argv);
                *(reinterpret_cast<double*>(argRec.data))= atof(argv[argument]);
                removeArgument(argument, argc, argv);
                break;

            case LONG:
                removeArgument(argument, argc, argv);
                *(reinterpret_cast<long*>(argRec.data)) =
                    atol(argv[argument]);
                removeArgument(argument, argc, argv);
                break;

            case STRING_LIST:  {
                removeArgument(argument, argc, argv);
                StringList* strList =
                    reinterpret_cast<StringList*>(argRec.data);
                // Add more entries until we hit the next parameter
                while ((argument < argc) && (argv[argument][0] != '-') &&
                       (argv[argument][1] != '-')) {
                    strList->push_back(argv[argument]);
                    removeArgument(argument, argc, argv);                    
                }
                break;
            }
                
            default:
                foundMatch = false;
                std::cerr << "Invalid arg type in arg array!" << std::endl;
                if (caxoe) {
                    exit(-1);
                }
            }
            if (foundMatch) {
                // we need to reprocess all of the args, to be safe...
                argument--;
            }
        }
    }
    if (caxoe) {
        checkRemainingArguments(argc, argv, true, false);
    }
}


void 
ArgParser::removeArgument(int removeIdx, int& argc, char* argv[]){
    memmove(argv + removeIdx, argv + removeIdx + 1,
            sizeof(char*) * (argc - removeIdx - 1));
    argv[--argc] = NULL;
}


bool 
ArgParser::checkRemainingArguments(int argc, char* argv[],
                                   bool caxoe, bool caxoa) {
    const std::string HelpStr("--help");
    for (int i = 1; (i < argc); i++) {
        if (HelpStr == argv[i]) {
            std::cout << *this << std::endl;
            if (caxoe) {
                exit(0);
            }
            return false;
        }
        
        if ((argv[i][0] == '-') || (caxoa)) {
            // Someone passed in an illegal argument!
            std::cerr << "Invalid argument(s) found. "
                      << "Offending arguments are:\n   ";
            for(int arg = 1; (arg < argc); arg++) {
                std::cerr << " " << argv[arg];
            }
            std::cerr << "\nSee help below for Valid arguments.\n\n";
            std::cerr << *this << std::endl;
            if (caxoe) {
                exit(-1);
            }
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, ArgParser& ap) {
    static bool helpFlag = true;
    static const ArgParser::ArgRecord HelpArgRec =
        {"--help", "Print this message", &helpFlag, ArgParser::BOOLEAN};
    const int numSpaces   = 3;
    const int indentation = 2;
    
    // Calculate the length of the longest argument.  Ensure that
    // maxLen is no smaller than 6 characters to handle the
    // omnipresent --help option
    int maxLen = 6;
    for (size_t i = 0; (i < ap.argRecords.size()); i++) {
        maxLen = std::max(maxLen, (int) ap.argRecords[i].command.size());
    }
    // Print the argument array
    const std::string IndentStr(indentation, ' ');
    for (size_t i = 0; (i < ap.argRecords.size()); i++) {
        const ArgParser::ArgRecord& argRec = ap.argRecords[i];
        if (argRec.type != ArgParser::MAIN_MESSAGE) {
            ap.printArg(os, IndentStr,
                        maxLen - argRec.command.size() + numSpaces,
                        argRec);
        } else {
            os << argRec.help << std::endl;
        }
    }
    // Print the omnipresent help option
    ap.printArg(os, IndentStr, maxLen + numSpaces - 6, HelpArgRec);
    return os;
}

std::string
ArgParser::getValue(const ArgRecord& argRec) {
    std::ostringstream retVal;
    switch(argRec.type) {
    case BOOLEAN:
        retVal << (*(reinterpret_cast<bool*>(argRec.data)) ? "true" : "false");
        break;
    case INTEGER:
        retVal << *(reinterpret_cast<int*>(argRec.data));
        break;
    case UNSIGNED_INT:
        retVal << *(reinterpret_cast<unsigned int*>(argRec.data));
        break;
    case FLOAT:
        retVal << *(reinterpret_cast<float*>(argRec.data));
        break;
    case DOUBLE:
        retVal << *(reinterpret_cast<double*>(argRec.data));
        break;
    case LONG:
        retVal << *(reinterpret_cast<long*>(argRec.data));
        break;                
    case STRING:
        retVal << *(reinterpret_cast<std::string*>(argRec.data));
        break;
    case STRING_LIST:
        // Currently we do nothing for string list.
        break;
    case INFO_MESSAGE:
        break;
    default:
        std::cerr << "Invalid arg type in arg array!" << std::endl;
    }
    // Return string representation from output stream
    return retVal.str();
}

void
ArgParser::printArg(std::ostream& os, const std::string& indentStr,
                    const int numSpaces, const ArgRecord& argRec,
                    const int maxWidth) {
    // Indent the proper amount
    os << indentStr;
    // here is the actual argument
    os << argRec.command;
    // print out the padding - leave numSpaces spaces between args and
    // help text...
    os << std::string(numSpaces, ' ');
    // here is the help string. Ensure it does not exceed max Width.
    const size_t spaceRemaining = maxWidth - indentStr.size() -
        argRec.command.size() - numSpaces - 1;
    const std::string defaultValue = getValue(argRec);
    std::string help = argRec.help + (!defaultValue.empty() ?
                                      " [" + defaultValue + "]" : "");
    while (!help.empty()) {
        os << help.substr(0, spaceRemaining);
        help = help.substr(std::min(spaceRemaining, help.size()));
        if (!help.empty()) {
            os << std::endl << std::string(maxWidth - spaceRemaining - 1, ' ');
        }
    }
    os << std::endl;
}

#endif