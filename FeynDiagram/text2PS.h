/*
 *------------------------------------------------------------------------
 * FeynDiagram
 *    This file is part of the FeynDiagram C++ library for creation
 *    of Feynman diagrams.
 *
 * Version $Revision: 1.2 $,  $Date: 2008/11/15 13:48:31 $
 *
 * by Bill Dimm,  bdimm@hotneuron.com
 *
 * Copyright 1993 Bill Dimm
 *    This code may be freely distributed.  No fee beyond a reasonable
 *    copying fee may be charged for this software.
 *    This code may not be included in commercial software without the
 *    author's permission.
 *------------------------------------------------------------------------
 */

#ifndef TEXT2PS_H
#  define TEXT2PS_H

// $Id: text2PS.h,v 1.2 2008/11/15 13:48:31 maitreda Exp $

#include <iostream>

extern void text2PS(std::ostream &out, char *textptr, double fontsize, char *fontname = "Times-Roman");

extern void text2PS(std::ostream &out, char *textptr, double horizfract, double vertfract, double fontsize, char *fontname = "Times-Roman");

#endif  // ifndef TEXT2PS_H
