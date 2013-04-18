ODFfibConvert
=============

To convert an ODF from nrrd format to fib format

###Usage
```
$ ./ODFfibConvert --help

USAGE: 
   ./ODFfibConvert  [--returnparameterfile <std::string>]
                    [--processinformationaddress <std::string>] [--xml]
                    [--echo] [--outputODF <std::string>] [--ODFfib
                    <std::string>] [--ODFitk <std::string>] [--]
                    [--version] [-h]
Where: 

   --returnparameterfile <std::string>
     Filename in which to write simple return parameters (int, float,
     int-vector, etc.) as opposed to bulk return parameters (image,
     geometry, transform, measurement, table).

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   --outputODF <std::string>
     Output ODF image in ITK or fib format

   --ODFfib <std::string>
     ODF image in fib format (.fib, .fib.gz)

   --ODFitk <std::string>
     ODF image in itk format (.nrrd, .nii, .nii.gz)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Description: Conversion ITK/fib. Give either an ODF image in itk format
   or fib format, if both are given the conversion ITK to fib will be
   done.

   Author(s): Adrien Kaiser

   Acknowledgements: Thank you everyone.
```

