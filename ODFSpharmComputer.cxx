
/* std classes */
#include <vector>
#include <fstream>
#include <iostream>

/* ITK classes */
#include <itkImage.h>
#include <itkImageBase.h> // for SizeType and SpacingType without having to declare an ImageType
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVectorImage.h>

/* Clapack for reconstruction of spharm */
extern "C" {
#include "f2c.h"
#include "clapack.h"
int sgels_(char *trans, integer *m, integer *n, integer *
           nrhs, real *a, integer *lda, real *b, integer *ldb, real *work, integer *lwork, integer *info);
}

#include "define.h" // types and constant definitions

void WriteOutITKODFImage( ODFImageType::Pointer ODFImage, std::string OutputODFImageFile ); // in ODFfibConvert.cxx

// spharm-pdm/Libraries/Shape/Algorithms/ParametricMeshToSPHARMSpatialObjectFilter.cxx: 999:  void ParametricMeshToSPHARMSpatialObjectFilter::ComputeCoeffs()
// called in                             ParametricMeshToSPHARMSpatialObjectFilter.cxx: 1233: void ParametricMeshToSPHARMSpatialObjectFilter::GenerateData()
// called at       spharm-pdm/Applications/ParaToSPHARMMeshCLP/ParaToSPHARMMeshCLP.cxx: 485:  spharmFilter->GenerateData();

  /////////////////////////////////////////
 //            LEGENDRE                 //
/////////////////////////////////////////

// Not modified
// Tabelliert P_m^m, P_{m+1}^m, ... P_l^m  in  p_ptr[0]... p_ptr[l-m].
// Das heisst p_ptr[ll-m] = P_ll^m,   fuer m <= ll <= l.
// Obiges bezieht sich auf den urspruenglich uebergebenen p_ptr.
void plgndr_row(const int l, const int m,
                const double x,
                const double somx2,             // somx2=sqrt((1.0-x)*(1.0+x));
                double *p_ptr)
{
  double pmm = 1.0;

  if( m > 0 )
    {
    double fact = 1.0;
    for( int i = 1; i <= m; i++ )
      {
      pmm *= -fact * somx2;
      fact += 2.0;
      }
    }
  *p_ptr++ = pmm;
  if( l > m )
    {
    double pmmp1 = x * (2 * m + 1) * pmm;
    *p_ptr++ = pmmp1;
    for( int ll = (m + 2); ll <= l; ll++ )
      {
      double pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
      *p_ptr++ = pll;
      pmm = pmmp1;
      pmmp1 = pll;
      }
    }
}

// Not modified
class ParametricMeshToSPHARMSpatialObjectFilterLegendre
{
  const double       SQRT2;
  double * * *       plm;

  double A(int l, int m)
  {
    int    i;
    double j;
  
    // (m >=0) && (l >= m)
    j = 1.0;
    for( i = l + m; i > l - m; i-- )
      {
      j *= (double)i;
      }
    return sqrt( (double) ( (l + l + 1) / j) / (4 * M_PI) );
  } // double A()

  double P(int l, int m, float z)
  {
    int    i;
    double f, p1o, p_1, p_2;

    if( m == 0 )
      {
      if( l == 0 )
        {
        return 1;
        }
      else
        {
        p_1 = z; p_2 = 1;
        for( i = 2; i <= l; i++ )
          {
          p1o = p_1;
          p_1 = ( (i + i - 1) * z * p_1 - (i - 1) * p_2) / double(i);
          p_2 = p1o;
          }
        return p_1;
        }
      }
    else
      {
      f = ( (plm[l])[m])[l - m];
      for( i = l - m - 1; i >= 0; i-- )
        {
        f = z * f + ( (plm[l])[m])[i];
        }
      return f;
      }
  } // double P()

  void ylm(int l, int m, float x, float y, float z, double & re, double & im)
  {
    int    i;
    double p,  p_ptr[256];

    switch( m )
      {
      case 0: re = 1.0; im = 0.0;
        break;
      case 1: re = x; im = y; // sincos(atan2(y, x), &im, &re);
        break;
      default:
        //
        re = (1.0 - z) * (1.0 + z);
        for( p = re, i = 1; i < m; i++ )
          {
          p *= re;
          }
        p = sqrt( (double) p);
        if( p < 1e-20 )
          {
          p = 0;
          }
        else
          {
          p = 1 / p;
          }
        im = sin(m * atan2(y, x) );
        re = cos(m * atan2(y, x) );
        re *= p; im *= p;
      }
    plgndr_row(l, m, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), p_ptr);
    p = p_ptr[l - m];
    i = l * (l + 1) / 2 + m;
    re *= p * a[i] * SQRT2; im *= p * a[i] * SQRT2;
  } // void ylm()

public:
  const int L;
  double *  a;
  void Ylm(int l, int m, float x, float y, float z, double & re, double & im)
  {
  
    if( m < 0 )
      {
      ylm(l, -m, x, y, z, re, im);
      if( m % 2 == 0 )
        {
        im = -im;
        }
      else
        {
        re = -re;
        }
      }
    else
      {
      ylm(l, m, x, y, z, re, im);
      }
    // cout << x << "," << y << "," << z << ": ";
    // cout << re << " " << im << ": ";
  } // void Ylm()

  void Ylm(int l, int m, float theta, float phi, double & re, double & im)
  {
    const int     MAX_L = 10000;
    static double p[MAX_L + 1];
    double        z;

    z = cos(theta);
    // sincos(m*phi, &im, &re);
    im = sin(m * phi);
    re = cos(m * phi);
    plgndr_row(l, m, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), p);
    re *= a[l * (l + 1) / 2 + m] * p[l - m] * SQRT2;
    im *= a[l * (l + 1) / 2 + m] * p[l - m] * SQRT2;
  } // void Ylm()

  ParametricMeshToSPHARMSpatialObjectFilterLegendre(const int l_in) :  SQRT2(sqrt( (double) 2) ),  L(l_in)
  {
    int    i, j, l, m, sgn;
    double n1, n2;

    // compute normalizing factors
    a = new double[(L + 1) + L * (L + 1) / 2];
    for( i = 0, l = 0; l <= L; l++ )
      {
      for( m = 0; m <= l; m++ )
        {
        a[i++] = A(l, m);
        }
      }
    // create associated Legendre polynomials
    // allocate array to array of polynomials
    plm = (double * * *) malloc( (L + 1) * sizeof(void *) ); // new double**[L+1];
    assert(plm != 0);
    for( l = 0; l <= L; l++ )
      {
      // allocate array of polynomials
      plm[l] = (double * *) malloc( (l + 1) * sizeof(void *) ); // new double*[l+1];
      assert(plm[l] != 0);
      for( m = 0; m <= l; m++ )
        {
        // allocate array for coeff's of polynomial
        (plm[l])[m] = (double *) malloc( (l - m + 2 + 1) * sizeof(double) );
        assert( (plm[l])[m] != 0);
        ( (plm[l])[m])[l - m + 1] = 0.0;
        ( (plm[l])[m])[l - m + 2] = 0.0;
        }
      }
    // P(0,0)
    ( (plm[0])[0])[0] = 1.0; ( (plm[0])[0])[1] = 0.0;  ( (plm[0])[0])[2] = 0.0;
    // P(1,0)
    ( (plm[1])[0])[0] = 0.0; ( (plm[1])[0])[1] = -1.0; // *(-1)^m !!!!!!!
    ( (plm[1])[0])[2] = 0.0; ( (plm[1])[0])[3] = 0.0;
    // P(1,1)
    ( (plm[1])[1])[0] = 1.0; ( (plm[1])[1])[1] = 0.0; ( (plm[1])[1])[2] = 0.0;
    for( l = 2; l <= L; l++ )
      {
      // compute coefficients of Legendre polynomials (recursive)
      n1 = (double)(l + l - 1); n2 = (double)(l - 1);
      // cout << "l" << l << ":";
      for( j = l; j > 0; j-- )
        {
        ( (plm[l])[0])[j] = (n1 * ( (plm[l - 1])[0])[j - 1] - n2 * ( (plm[l - 2])[0])[j])
          / (double)l;
        // cout << " " << ((plm[l])[0])[j];
        }
      ( (plm[l])[0])[0] = -n2 * ( (plm[l - 2])[0])[0] / (double)l;
      // cout << " " << ((plm[l])[0])[0];
      sgn = 1;
      for( m = 1; m <= l; m++ )
        {
        // differentiate Legendre polynomials m times
        // cout << " m" << m << ":";
        sgn *= -1;
        for( j = l - m; j >= 0; j-- )
          {
            ( (plm[l])[m])[j] = sgn * (j + 1) * ( (plm[l])[m - 1])[j + 1];
          // cout << " " << ((plm[l])[m])[j];
          }
        }
      // cout << endl;
      }
  } // ParametricMeshToSPHARMSpatialObjectFilterLegendre()

  ~ParametricMeshToSPHARMSpatialObjectFilterLegendre()
  {
    int l, m;
    for( l = 0; l <= L; l++ )
      {
      for( m = 0; m <= l; m++ )
        {
        free( (plm[l])[m]);
        }
      free(plm[l]);
      }
    free(plm);
    delete [] a;
    // cout << "destructor called!"<<endl;
  } // ~ParametricMeshToSPHARMSpatialObjectFilterLegendre()
};


  /////////////////////////////////////////
 //            FUNCTIONS                //
/////////////////////////////////////////

// Modified
int Get_BaseVal( float * & A, int & nvert, int & degree )
{
  // Added:

  // get odf vertices from file
  std::ifstream InfileStream (DATFILE , std::ios::in); // open in reading

  std::string line;
  double xVert, yVert, zVert;

  std::getline( InfileStream, line ); // void read to avoid reading the first line as point

  nvert = NBDIRS; // 321
  double * vert = new double[nvert * 3];
  int i=0;
  while( std::getline( InfileStream, line ) && i<nvert ) 
  {
    std::stringstream linestr(line);
    linestr >> xVert >> yVert >> zVert;

    vert[3 * i + 0] = xVert;
    vert[3 * i + 1] = yVert;
    vert[3 * i + 2] = zVert;

    i++;
  }

  InfileStream.close();

  // Original Get_BaseVal():

  double re, im, z;

  unsigned int m_Degree = SPHARMDEGREE; // 4
  ParametricMeshToSPHARMSpatialObjectFilterLegendre * m_leg;

  if( m_leg )
    {
    delete m_leg;
    }
  m_leg = new ParametricMeshToSPHARMSpatialObjectFilterLegendre(m_Degree);

  double * plm = new double[(m_Degree + 1) * (m_Degree + 1) * 3]; // 25 * 3

//  degree = m_leg->L + 1;
//  degree *= degree; // fuer re und im
  degree = NBCOEFFS; // 15
  A = new float[degree * nvert];
  // for all vertices nvert = NBDIRS = 321
  for( int ind = 0, i = 0; i < nvert; i++, ind += 3 )
  {
    z = vert[ind + 2];

    // void plgndr_row(const int l, const int m, const double x, const double somx2, double *p_ptr) // somx2=sqrt((1.0-x)*(1.0+x));
    // Tabelliert P_m^m, P_{m+1}^m, ... P_l^m  in  p_ptr[0]... p_ptr[l-m].
    // Das heisst p_ptr[ll-m] = P_ll^m,   fuer m <= ll <= l.

    // Ylm: Remove odd coefficients: l=1, l=3: keep l = 0, 2, 4
    unsigned int Lindextable[3] = { 0, 1, 6 };

    // Yl0
    plgndr_row(m_leg->L, 0, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), plm);
    // m_leg->L = m_Degree = SPHARMDEGREE = 4
//    for( int l = 0; l <= m_leg->L; l++ )
    for( int l = 0; l <= m_leg->L; l+=2 ) // l = 0, 2, 4
    {
//      A[i + l * l * nvert] = m_leg->a[l * (l + 1) / 2] * plm[l];      // assign Y_l^0
      A[i + Lindextable[l/2] * nvert] = m_leg->a[l * (l + 1) / 2] * plm[l];      // assign Y_l^0
    }

    // Yl1 Yl-1
    plgndr_row(m_leg->L, 1, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), plm);
    im = sin(atan2(vert[ind + 1], vert[ind]) );
    re = cos(atan2(vert[ind + 1], vert[ind]) );
    // m_leg->L = m_Degree = SPHARMDEGREE = 4
//    for( int l = 1; l <= m_leg->L; l++ )
    for( int l = 2; l <= m_leg->L; l+=2 ) // l = 2, 4
    {
//      int j = l * l + 1;
      int j = Lindextable[l/2] + 1;
      A[i + j * nvert] = m_leg->a[l * (l + 1) / 2 + 1] * plm[l - 1] * re;
      A[i + j * nvert + nvert] = m_leg->a[l * (l + 1) / 2 + 1] * plm[l - 1] * im;
    }

    // Yl2 Yl-2 .. YlL Yl-L
    // m_leg->L = m_Degree = SPHARMDEGREE = 4
    for( int m = 2; m <= m_leg->L; m++ )
    {
      im = sin(m * atan2(vert[ind + 1], vert[ind]) );
      re = cos(m * atan2(vert[ind + 1], vert[ind]) );
      // double p= (1.0-z)*(1.0+z);
      // p= sqrt(pow(p, m));
      // if (p<1e-20) p= 0; else p= 1/p;
      // re*= p; im*= p;
      plgndr_row(m_leg->L, m, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), plm);
      // m_leg->L = m_Degree = SPHARMDEGREE = 4
      for( int l = m; l <= m_leg->L; l++ )
      {
        if( l%2 == 0 ) // if l is even: OK
        {
//        int j = l * l + m + m - 1;
          int j = Lindextable[l/2] + m + m - 1;
          A[i + j * nvert] = m_leg->a[l * (l + 1) / 2 + m] * plm[l - m] * re;
          A[i + j * nvert + nvert] = m_leg->a[l * (l + 1) / 2 + m] * plm[l - m] * im;
        } // if l is even
      }
    }
  }

  delete plm;
  delete vert;

  return 0;
}

// Homemade
ODFType *OrderCoeffs (float* obj)
{
/*
obj: 25 coeffs
0    1   2   3   4    5   6    7   8    9   10  11   12  13   14  15  16  17  18   19   20   21  22  23  24
C00 C10 C11 C1-1 C20 C21 C2-1 C22 C2-2 C30 C31 C3-1 C32 C3-2 C33 C3-3 C40 C41 C4-1 C42 C4-2 C43 C4-3 C44 C4-4
-> spharm-pdm/Libraries/Shape/Numerics/SphericalHarmonicsPolynomial.txx:56:Evaluate()

Updated:
obj: 15 coeffs
0    1   2   3   4    5   6    7   8    9   10  11   12  13   14
C00 C20 C21 C2-1 C22 C2-2 C40 C41 C4-1 C42 C4-2 C43 C4-3 C44 C4-4
    ||
    VV
coeffsOrdered: 15 coeffs (no odd orders/degrees) (same order as in cost)
0     1    2   3  4    5   6   7     8    9   10   11  12  13  14
C00 C2-2 C2-1 C20 C21 C22 C4-4 C4-3 C4-2 C4-1 C40 C41 C42 C43 C44
-> COST/ODFReconstructor.cxx:136:Y()
                            :174:GetLM()
*/

  ODFType *coeffsOrdered = new ODFType[NBCOEFFS];

  coeffsOrdered[  0 ] = obj[  0 ]; //  0 ];
  coeffsOrdered[  1 ] = obj[  5 ]; //  8 ];
  coeffsOrdered[  2 ] = obj[  3 ]; //  6 ];
  coeffsOrdered[  3 ] = obj[  1 ]; //  4 ];
  coeffsOrdered[  4 ] = obj[  2 ]; //  5 ];
  coeffsOrdered[  5 ] = obj[  4 ]; //  7 ];
  coeffsOrdered[  6 ] = obj[ 14 ]; // 24 ];
  coeffsOrdered[  7 ] = obj[ 12 ]; // 22 ];
  coeffsOrdered[  8 ] = obj[ 10 ]; // 20 ];
  coeffsOrdered[  9 ] = obj[  8 ]; // 18 ];
  coeffsOrdered[ 10 ] = obj[  6 ]; // 16 ];
  coeffsOrdered[ 11 ] = obj[  7 ]; // 17 ];
  coeffsOrdered[ 12 ] = obj[  9 ]; // 19 ];
  coeffsOrdered[ 13 ] = obj[ 11 ]; // 21 ];
  coeffsOrdered[ 14 ] = obj[ 13 ]; // 23 ];

  return coeffsOrdered;
}

// Modified
// Compute all 25 coeffs for 1 voxel only
ODFType *ComputeCoeffs( ODFType* ODFSampledImageArray, unsigned long VoxelIndex, float * A, int m_int, int n_int )
{
/* Solving system A.X = B
A = basis matrix depending on the vertices only (computed once at the beginning)  321 rows x 25 cols
X = 25 shparm coeffs 25 rows x 1 col
B = ODF values only  321 rows x 1col
*/

  integer numRH, m, n;
  m = m_int;
  n = n_int;
  numRH = 1;
  float * obj = new float[m * numRH]; // will contain data and then overwritten with results
  // Get odf data in obj array for equation solving
  for( int i = 0; i < m; i++ ) // m will be NBDIRS
  {
    obj[i] = ODFSampledImageArray[ VoxelIndex*NBDIRS + i ];
  }

  char    trans[20] = "N"; // Use B without transposing it
  int     fac = n * m;
  integer workSize = fac * 2;
  integer info; // return code
  float * work = new float[workSize]; // workspace

  integer lda = m;
  integer ldb = m;

  sgels_(trans, &m, &n, &numRH, A, &lda, obj, &ldb, work, &workSize, &info); // 'integer' type needed for sgels_()
/* From LAPACK: Solve G(eneral) Equations Linear System _  // http://www.netlib.org/lapack/single/sgels.f
   Solve A.X = B
   trans   = 'N' or 'T': Use B as is or transposed
   m       = nb rows in A
   n       = nb columns in A
   numRH   = nb of Right Hand sides = nb of columns in B (only 1 here)
   A       = real array of dim (lda, n) // A is modfied by sgels_
   lda     = max(1,m)
   obj = B = real array of dim (ldb,numRH)
   On exit, B is overwritten by the solution vectors, stored columnwise:
   rows 1 to n of B contain the least squares solution vectors //  ==> obj contains the flattened coefs
   ldb     = max(1,m,n)
   work    = workspace
   workSize= size of work
   info    = return code (0, >0 or <0)
*/

  delete work; // new done in ComputeCoeffs()
  work = NULL;

  ODFType *coeffsOrdered = OrderCoeffs(obj);

  delete [] obj;
  obj = NULL;
  
  return coeffsOrdered;
}

// Used with forward declaration using function prototype in ODFfibConvert.cxx
// Compute coeffs for all voxels and write out image
void ComputeSpharmCoeffs( ODFImageType::Pointer ODFSampledImage,
                          std::string OutputODFCoeffs,
                          itk::ImageBase< 3 >::SizeType size,
                          itk::ImageBase< 3 >::SpacingType spacing)
{
  std::cout<< "Status : Computing Spherical Harmonics Coefficients from ODF samples"<< std::endl;

  //// Set Image properties
  ODFImageType::Pointer NewODFCoeffsImage = ODFImageType::New();
  NewODFCoeffsImage->SetSpacing ( spacing ) ;
  NewODFCoeffsImage->SetVectorLength ( NBCOEFFS ) ;
  // origin
  ODFImageType::PointType origin ;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  NewODFCoeffsImage->SetOrigin ( origin ) ;
  // region (size)
  ODFImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  ODFImageType::RegionType region;
  region.SetSize ( size ) ;
  region.SetIndex ( start ) ;
  NewODFCoeffsImage->SetRegions ( region ) ;
  // Allocate image
  NewODFCoeffsImage->Allocate () ;
  NewODFCoeffsImage->FillBuffer ( 0 ) ;
  // Buffer Array
  ODFType *NewODFCoeffsImageArray;
  NewODFCoeffsImageArray = NewODFCoeffsImage->GetPixelContainer()->GetBufferPointer() ;

  // Compute Sherical Harmonics coefficients

  // Compute basis matrix A: TO BE DONE ONLY ONCE FOR ALL VOXELS (depends ONLY on the odf vertices)
  float * A;
  int NbVertices;
  int NbCoeffs;
  // Compute basis matrix A: depends on vertices ONLY
  Get_BaseVal( A, NbVertices, NbCoeffs ); // will set A, NbVertices, NbCoeffs (reference &)

  ODFType* ODFSampledImageArray = ODFSampledImage->GetPixelContainer()->GetBufferPointer();

  // for all voxels in the image
  int percent=0;
  for( unsigned long VoxelIndex = 0 ; VoxelIndex < size[0]*size[1]*size[2] ; VoxelIndex++ )
  {
    // display every 10%
    if( VoxelIndex % (size[0]*size[1]*size[2]/10) == 0 ) // does not work if nbVoxels < 10 (OK for 3D images)
    {
      std::cout<< percent << "%" <<std::endl;
      percent += 10;
    }

    // create a new copy of A as it has been modified by sgels_()
    float * Acopy = new float[ NbVertices * NbCoeffs ];
    for( unsigned int Aindex=0 ; Aindex < NbVertices*NbCoeffs ; Aindex++ )
    {
      Acopy[ Aindex ] = A[ Aindex ];
    }

    ODFType *SpharmCoeffs = ComputeCoeffs( ODFSampledImageArray, VoxelIndex, Acopy, NbVertices, NbCoeffs ); // first 25 values of array are results = sph. harm. coeffs

    delete [] Acopy;
    Acopy = NULL;

    for( unsigned int coeff=0; coeff < NBCOEFFS ; ++coeff ) // for all coeffs
    {
      NewODFCoeffsImageArray[ VoxelIndex*NBCOEFFS + coeff ] = SpharmCoeffs[ coeff ];
    } // for all coeffs

    delete [] SpharmCoeffs; // new done in OrderCoeffs(): ODFType *coeffsOrdered = new ODFType[]; ... return coeffsOrdered;
    SpharmCoeffs = NULL;

  } // for all voxels in the image
  std::cout<< "100%" <<std::endl;

  delete A; // new done in Get_BaseVal(), called once so 1 delete
  A = NULL;

  // write out coeffs image
  WriteOutITKODFImage( NewODFCoeffsImage, OutputODFCoeffs );

  return ;
}

