ODFfibConvert
=============

To convert an ODF from nrrd format to fib format

###Usage
```
$ ./ODFfibConvert --help

USAGE: 

   ./ODFfibConvert  [--returnparameterfile <std::string>]
                    [--processinformationaddress <std::string>] [--xml]
                    [--echo] [--outputODF <std::string>] [--fa
                    <std::string>] [--mask <std::string>] [--ODFfib
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

   --fa <std::string>
     FA image in itk format (.nrrd, .nii, .nii.gz)

   --mask <std::string>
     Mask image to apply to the ODF in itk format (.nrrd, .nii, .nii.gz)

   --ODFfib <std::string>
     ODF image in fib format (.fib, .fib.gz)

   --ODFitk <std::string>
     ODF image in itk format (.nrrd, .nii, .nii.gz). A mask and a FA image
     are required to convert ITK image to fib file.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Description: Conversion ITK/fib. Give either an ODF image in itk format
   or fib format, if both are given the conversion ITK to fib will be done.
   The ODF itk image must be a vector image containing the 15 spherical
   harmonics coefficients for each voxel.

   Author(s): Adrien Kaiser

   Acknowledgements: Thank you everyone.
```

