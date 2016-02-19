/*
 * Global_constants.h
 *
 *  Created on: Jun 3, 2015
 *      Author: drasal
 */
#ifndef INCLUDE_GLOBAL_CONSTANTS_H_
#define INCLUDE_GLOBAL_CONSTANTS_H_

#include <string>
#include <vector>
#include <math.h>

namespace analyzer {

  static const std::string delphesTreeName  = "DelphesSim";

  static const int         proton           = 2212;
  static const int         gamma            =   22;
  static const int         b                =    5;
  static const int         bbar             =   -5;
  static const int         t                =    6;
  static const int         tbar             =   -6;
  static const int         higgs            =   25;
  static const int         Z0               =   23;
  static const int         WPls             =   24;
  static const int         WMin             =  -24;
  static const int         electron         =   11;
  static const int         positron         =  -11;
  static const int         muon             =   13;
  static const int         antimuon         =  -13;
  static const int         tau              =   15;
  static const int         antitau          =  -15;
  static const int         undefined        =    0;
}

#endif /* INCLUDE_GLOBAL_CONSTANTS_H_ */
