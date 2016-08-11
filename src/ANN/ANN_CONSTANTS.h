//----------------------------------------------------------------------
// File: ANN_CONSTANTS.h
// Programmer:		Constantino A. Garcia
// Description:		Basic include file with ANN constants
// Last modified:	19/07/2015
#ifndef ANN_CONSTANTS_H
#define ANN_CONSTANTS_H
 
// define your own namespace to hold constants
namespace ann_constants{
   /* constants used in kd_dump.cpp */
  const int  	ANN_STRING_LEN		= 500;	// maximum string length
  const double	ANN_EPSILON			= 1E-5; // small number for float comparison
  /* constants used in kd_split.cpp */
  const double ANN_ERR = 0.001;				// a small value
  const double ANN_FS_ASPECT_RATIO = 3.0;		// maximum allowed aspect ratio
	                    									// in fair split. Must be >= 2.

}

#endif
