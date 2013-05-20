/*
./ODFfibConvert --ODFfib /home/akaiser/Networking/NFG/nfg_1.1.1/generated_collections/130325112527/DWI/dwi-00.src.gz.odf8.f5rec.qbi.sh8.0.006.fib.gz --outputSampledODF ./ODFsampled.nrrd
./ODFfibConvert --ODFitk /home/akaiser/Networking/NFG/nfg_1.1.1/generated_collections/130325112527/DWI/odf.nrrd --mask /home/akaiser/Networking/NFG/nfg_1.1.1/generated_collections/130325112527/DWI/maskDSI.nii.gz --fa /home/akaiser/Networking/NFG/nfg_1.1.1/generated_collections/130325112527/DWI/fa.nrrd --outputSampledODF ./ODFsampled.nrrd
*/

/* std classes */
#include <vector>
#include <fstream>
#include <iostream>

#include <zlib.h>

/* ITK classes */
#include <itkImage.h>
#include <itkImageBase.h> // for SizeType and SpacingType without having to declare an ImageType
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVectorImage.h>
#include <itksys/SystemTools.hxx>

/* DSIstudio header */
#include <libs/mat_file.hpp>
#include <libs/utility/prog_interface.cpp>

/* Homemade classes and headers */
#include "ODFfibConvertCLP.h" //generated when ccmake
#include "ODFReconstructor.h" // Comes from COST, modified to have ODFType = float instead of double and comment all std::cout
#include "define.h" // types and constant definitions

// ODFSpharmComputer.cxx
void ComputeSpharmCoeffs( ODFImageType::Pointer ODFSampledImage,
                          std::string OutputODFCoeffs,
                          itk::ImageBase< 3 >::SizeType size,
                          itk::ImageBase< 3 >::SpacingType spacing);

  /////////////////////////////////////////
 //              COMMON                 //
/////////////////////////////////////////

void WriteOutITKODFImage( ODFImageType::Pointer ODFImage, std::string OutputODFImageFile )
{
  std::cout<<"Status : Writing out ITK ODF image: " << OutputODFImageFile << " ... ";

  typedef itk::ImageFileWriter < ODFImageType > ODFWriterType ;
  ODFWriterType::Pointer Writer = ODFWriterType::New() ;
  Writer->SetFileName ( OutputODFImageFile ); 
  Writer->SetInput( ODFImage );
  Writer->SetUseCompression(true);

  try
  {
    Writer->Update();
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Problem writing the image: " << OutputODFImageFile << std::endl ;
    std::cerr << excp << std::endl ;
  }

  std::cout<<"DONE"<<std::endl;

  return;
}

  /////////////////////////////////////////
 //            ITK -> FIB               //
/////////////////////////////////////////

template< class T > T* GetArrayFromFile( std::string ArrayFile, bool SkipFirstLine ) // SkipFirstLine needed for vertices file ( = datFile -> first line = nb of vertices)
{
  std::ifstream InfileStream (ArrayFile.c_str() , std::ios::in); // open in reading
  if(! InfileStream) // error while opening
  {
    return NULL;
  }

  // Get data in a 2 dim vector
  std::vector< std::vector< T > > Vertices;
  Vertices.resize(3);

  std::string line;
  T x, y, z;

  if( SkipFirstLine ) // void read to avoid reading the first line as triple
  {
    std::getline( InfileStream, line );
  }

  while( std::getline( InfileStream, line ) ) 
  {
    std::stringstream linestr(line);
    linestr >> x >> y >> z;

    Vertices[0].push_back(x);
    Vertices[1].push_back(y);
    Vertices[2].push_back(z);
  }

  InfileStream.close();

  // Write data in a 1 dim array
  T* data = new T [ 3 * Vertices[0].size() ] ; // delete after function call in writeFib()

  for( int i = 0 ; i < 3 ; i++ )
  {
    for( int j = 0 ; j < Vertices[0].size() ; j++ )
    {
      data[ i*Vertices[0].size() + j ] = Vertices[i][j];
    }
  }

  return data;
}

ODFImageType::Pointer loadITKODF(std::string ODFITKImageFile, std::string outputSampledODF ) // returns the sampled ODF image
{
  typedef itk::ImageFileReader< ODFImageType >  ODFReaderType;
  ODFReaderType::Pointer ODFReader       = ODFReaderType::New();
  ODFImageType::Pointer  ODFCoeffsImage  = ODFImageType::New();

  std::cout<<"Status : Loading ODF ITK image: "<< ODFITKImageFile <<" ... ";
  ODFReader->SetFileName ( ODFITKImageFile ) ;
  try
  {
    ODFReader->Update() ; 
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cout << "FAIL" << std::endl << "Error  : Problem reading the input ODF file: " << ODFITKImageFile <<  std::endl;
    std::cout << excp << std::endl;
    return NULL;
  }
  ODFCoeffsImage = ODFReader->GetOutput();
  std::cout<<"DONE"<<std::endl;

  if( ODFCoeffsImage->GetVectorLength() == 15 ) // ODF is given with coeffs
  {
    // Reconstruct ODF: convert sph. harm. coeffs to ODF value by vertex by sampling coeffs on vertices
    std::cout<<"Status : Reconstructing ODF (can take some time)"<<std::endl;
    ODFReconstructor ODFreconstructor ( ODFCoeffsImage, NBDIRS, 15, DATFILE ); // ODFReconstructor(ODFImageType::Pointer Input_image, unsigned long numberOfSamplesOnSphere, unsigned long numberOfSpharm, std::string filename )
    ODFImageType::Pointer ODFSampledImage = ODFreconstructor.ReconstructODFImage(); // ODFSampledImage is a vector image that contains sampled ODF for each vertex on the sphere
    // !! Nb of components is NOT doubled here (see ODFReconstructor.cxx:326)

    if( outputSampledODF != "" )
    {
      WriteOutITKODFImage( ODFSampledImage, outputSampledODF );
    }

//  ComputeSpharmCoeffs( ODFSampledImage, "/home/akaiser/work/Projects/ODFfibConvert-build/ODFcoeffsFromReconstructed2.nrrd", ODFSampledImage->GetLargestPossibleRegion().GetSize(), ODFSampledImage->GetSpacing() );

    return ODFSampledImage;
  }
  else if( ODFCoeffsImage->GetVectorLength() == 321 ) // ODF is given sampled with the right nb of vertices
  {
    return ODFCoeffsImage; // ODFCoeffsImage is actually sampled
  }
  else
  {
    std::cout <<"Error  : Input ODF image must have either 15 (spharm coeffs) or 321 (ODF samples) components."<<std::endl;
    return NULL;
  }

} // loadITKODF()

template< class T > typename itk::Image< T, 3 >::Pointer loadScalarITKImage( std::string ScalarImageFilename )
{
  typedef itk::Image < T, 3 > ImageType ;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typename ReaderType::Pointer Reader = ReaderType::New();
  typename ImageType::Pointer  Image  = ImageType::New();

  std::cout<<"Status : Loading scalar ITK image: "<< ScalarImageFilename <<" ... ";

  Reader->SetFileName ( ScalarImageFilename ) ;
  try
  {
    Reader->Update() ; 
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Error  : Problem reading the input scalar image: " << ScalarImageFilename <<  std::endl;
    std::cerr << excp << std::endl;
    return NULL;
  }
  Image = Reader->GetOutput();

  std::cout<<"DONE"<<std::endl;

  return Image;
}

unsigned long MaskFAImage( double* FAArray, short* MaskArray, unsigned long ArraySize )
{
  unsigned long NbVoxelsInMask = 0;

  for( unsigned long index = 0 ; index < ArraySize ; index++ )
  {
    if( MaskArray[index] == 0 ) // if outside of mask
    {
      FAArray[index] = 0;
    }
    else
    {
      NbVoxelsInMask++;
    }
  }
  return NbVoxelsInMask;
}

unsigned int WriteODFBlock (unsigned int BlockSize,
                            unsigned int ODFBlockIndex,
                            unsigned int Voxelindex,
                            double* FAArray,
                            short* MaskArray,
                            ODFType* ODFArray,
                            MatFile* mat_writer) // returns updated Voxelindex
{
/*
ODFArray (ITK reconstructed image):
  321   321   321   321
|-----|-----|-----|-----|... x 1000000 voxels

ODFFibArray (Fib image):
  321   321   321   321
|-----|-----|-----|-----|... x 20000 voxels   x  N blocks

dir        E [0 ; 320]    -> dir
BlockVoxel E [0 ; 19999]  -> ODFBlockVoxelIndex
voxel      E [0 ; 999999] -> Voxelindex

==> ODFFibArray[ BlockVoxel*321 + dir ] = ODFArray[ voxel*642 + dir ] ;
*/
  std::vector< ODFType > ODFFibVector;
  ODFFibVector.resize( BlockSize * NBDIRS );

  // for all voxels in the block
  unsigned int ODFBlockVoxelIndex=0;
  while ( ODFBlockVoxelIndex < BlockSize )
  {
    if( FAArray[ Voxelindex ] != 0.0 ||  MaskArray[ Voxelindex ] != 0.0 ) // if in mask // Inside of mask == at least one of fa0 or index0 matrices show != 0
    {
      for ( unsigned int dir=0; dir < NBDIRS ; ++dir ) // for all directions
      {
        ODFFibVector[ ODFBlockVoxelIndex*NBDIRS + dir ] = ODFArray[ Voxelindex*NBDIRS + dir ];
      } // for all directions

      ODFBlockVoxelIndex++;

    } // if in mask

    Voxelindex++;

  } // for all voxels in the block

  // Write out block
  std::ostringstream ODFBlockIndexout;
  ODFBlockIndexout << ODFBlockIndex;
  std::string ODFMatrixName = "odf" + ODFBlockIndexout.str();

  mat_writer->add_matrix(ODFMatrixName.c_str(), &*ODFFibVector.begin(), NBDIRS, BlockSize);

  return Voxelindex;
} // unsigned int WriteODFBlock()


// Master function

bool writeFib (ODFImageType::Pointer ODFSampledImage,
               std::string outputODFfib,
               std::string maskImageFilename,
               std::string FAImageFilename)
{

  MatFile mat_writer(outputODFfib.c_str());

// add_matrix creates the matrix and write it to file
// DSIstudio/libs/mat_file.hpp:414 : template<typename Type> void add_matrix(const char* name,const Type* data_ptr,unsigned int rows,unsigned int cols)

  //Size
  ODFImageType::SizeType size  = ODFSampledImage->GetLargestPossibleRegion().GetSize();
  short sizeArray[3] = {size[0], size[1], size[2]}; // SizeValueType = unsigned long  -> short
  mat_writer.add_matrix("dimension", sizeArray, 1, 3);

  // Spacing
  ODFImageType::SpacingType spacing = ODFSampledImage->GetSpacing();
  float spacingArray[3] = {spacing[0], spacing[1], spacing[2]}; // SpacingValueType = double  -> float
  mat_writer.add_matrix("voxel_size", spacingArray, 1, 3);

  // odf_vertices
  double* VerticesArray = GetArrayFromFile< double >(DATFILE, true);
  mat_writer.add_matrix("odf_vertices", VerticesArray, 3, NBDIRS*2);
  delete []VerticesArray;

  // odf_faces
  short* FacesArray = GetArrayFromFile< short >(FACESFILE, false);
  mat_writer.add_matrix("odf_faces", FacesArray, 3, NBFACES);
  delete []FacesArray;

// TODO: Compare size of images (ODF, mask and fa) to check they are the same

  // Mask as index0 (for mask reconstruction)
  itk::Image< short, 3 >::Pointer  MaskImage = loadScalarITKImage< short >( maskImageFilename );
  short* MaskArray = MaskImage->GetPixelContainer()->GetBufferPointer();
  mat_writer.add_matrix( "index0", MaskArray, 1, size[0]*size[1]*size[2] ); // DSI studio seg fault when loading if no index0 (needed for reconstruction of mask)

  // FA
  // !! declaration of itk::Image< T, 3 >::Pointer needed otherwise empty arrays
  itk::Image< double, 3 >::Pointer FAImage   = loadScalarITKImage< double >( FAImageFilename );
  double* FAArray  = FAImage->GetPixelContainer()->GetBufferPointer();
  unsigned long NbVoxelsInMask = MaskFAImage( FAArray, MaskArray, size[0]*size[1]*size[2] ); // will modify FAArray
  mat_writer.add_matrix( "fa0", FAArray, 1, size[0]*size[1]*size[2] );

  // Read ODF array and fill vector with matrices to be written out
  std::cout<<"Status : Converting ODF data" <<std::endl;
  ODFType* ODFArray = ODFSampledImage->GetPixelContainer()->GetBufferPointer();
  unsigned int NbBlocks = ( NbVoxelsInMask - NbVoxelsInMask%ODFBLOCKSIZE ) / ODFBLOCKSIZE ; // blocks of 20000 without the eventual last block
  unsigned int NbVoxelsLastBlock = NbVoxelsInMask%ODFBLOCKSIZE ;
  unsigned long Voxelindex = 0 ;

  for( unsigned int ODFBlockIndex=0; ODFBlockIndex < NbBlocks; ++ODFBlockIndex ) // for all 20000 voxels blocks
  {
    Voxelindex = WriteODFBlock(ODFBLOCKSIZE, ODFBlockIndex, Voxelindex, FAArray, MaskArray, ODFArray, &mat_writer);

    std::cout<< ODFBlockIndex+1 << "/" << NbBlocks + (int)(NbVoxelsLastBlock!=0) <<" blocks written" <<std::endl;
  }

  if( NbVoxelsLastBlock != 0 ) // if voxels in the last block
  {
    Voxelindex = WriteODFBlock(NbVoxelsLastBlock, NbBlocks, Voxelindex, FAArray, MaskArray, ODFArray, &mat_writer);
    NbBlocks++; // for display
    std::cout<< NbBlocks << "/" << NbBlocks <<" blocks written" <<std::endl;
  }

  std::cout<<"Status : ODF converted in "<< NbBlocks <<" blocks" <<std::endl;

  std::cout<<"Status : Done writing out fib file: " << outputODFfib << std::endl;

  return true;
}

  /////////////////////////////////////////
 //            FIB -> ITK               //
/////////////////////////////////////////

template<class T > void WriteITKScalarImage( T* data_ptr,
                                             std::string OutputITKScalar,
                                             typename itk::ImageBase< 3 >::SizeType size,
                                             typename itk::ImageBase< 3 >::SpacingType spacing )
{
  std::cout<<"Status : Writing out ITK scalar image: " << OutputITKScalar << " ... ";

  // ITK types and definitions
  typedef itk::Image <T, 3>  ImageType ;
  typename ImageType::Pointer NewImage = ImageType::New() ;
  NewImage->SetSpacing ( spacing ) ;

  //// Set Image properties

  // origin
  typename ImageType::PointType origin ;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  NewImage->SetOrigin ( origin ) ;

  // region (size)
  typename ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  typename ImageType::RegionType region;
  region.SetSize ( size ) ;
  region.SetIndex ( start ) ;
  NewImage->SetRegions ( region ) ;

  // Allocate image
  NewImage->Allocate () ;
  NewImage->FillBuffer ( 0 ) ;

  // Buffer Array
  T *NewImageArray;
  NewImageArray = NewImage->GetPixelContainer()->GetBufferPointer() ;

  // Fill buffer
  for( unsigned long i = 0 ; i < size[0]*size[1]*size[2] ; i++ )
  {
    NewImageArray[ i ] = data_ptr[ i ];
  }

  // Write out image
  typedef itk::ImageFileWriter < ImageType > WriterType ;
  typename WriterType::Pointer Writer = WriterType::New() ;
  Writer->SetFileName ( OutputITKScalar ); 
  Writer->SetInput( NewImage );
  Writer->SetUseCompression(true);

  try
  {
    Writer->Update() ;
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Problem writing the image: " << OutputITKScalar << std::endl ;
    std::cerr << excp << std::endl ;
  }

  std::cout<<"DONE"<<std::endl;

  return;
} // void WriteITKScalarImage()

template<class T > void GetVertices( T* data_ptr, int nbVertices, std::string VerticesOutputFile ) // matrix dim = (3, nbVertices)
{
  std::cout<<"Status : Writing file: " << VerticesOutputFile <<std::endl;

  std::ofstream VerticesFileStream (VerticesOutputFile.c_str() , std::ios::out | std::ios::trunc); // opening in writing with erasing the open file
  if(! VerticesFileStream)
  {
    std::cout<<"Error : Creating file: "<< VerticesOutputFile <<std::endl;
    return;
  }

  for( int j = 0 ; j < nbVertices ; j++ ) // 1 vertex by line
  {
    VerticesFileStream << data_ptr[ j ] << " " << data_ptr[ nbVertices + j ] << " " << data_ptr[ 2*nbVertices + j ] << std::endl;
  }

  VerticesFileStream.close();

  return;
}

ODFImageType::Pointer AssembleODFBlocks( std::vector< ODFType* > ODFBlocks,
                                         std::vector< unsigned int > BlockSizes,
                                         std::string outputSampledODF,
                                         itk::ImageBase< 3 >::SizeType size,
                                         itk::ImageBase< 3 >::SpacingType spacing,
                                         double* fa,
                                         short* MaskArray )
{
  std::cout<< "Status : Assembling ODF blocks to create ODF samples image"<< std::endl;

  // ITK types and definitions
  ODFImageType::Pointer NewODFSampledImage = ODFImageType::New() ;

  //// Set Image properties
  NewODFSampledImage->SetSpacing ( spacing ) ;
  NewODFSampledImage->SetVectorLength ( NBDIRS ) ;

  // origin
  ODFImageType::PointType origin ;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  NewODFSampledImage->SetOrigin ( origin ) ;

  // region (size)
  ODFImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  ODFImageType::RegionType region;
  region.SetSize ( size ) ;
  region.SetIndex ( start ) ;
  NewODFSampledImage->SetRegions ( region ) ;

  // Allocate image
  NewODFSampledImage->Allocate () ;
  NewODFSampledImage->FillBuffer ( 0 ) ;

  // Buffer Array
  ODFType *NewODFSampledImageArray;
  NewODFSampledImageArray = NewODFSampledImage->GetPixelContainer()->GetBufferPointer() ;

  // Fill buffer by assembling all ODF arrays
/*
ODFArray (ITK image):
  321   321   321   321
|-----|-----|-----|-----|... x 1000000 voxels

ODFFibArray (Fib image):
  321   321   321   321
|-----|-----|-----|-----|... x 20000 voxels   x  N blocks

dir        E [0 ; 320]    -> dir
BlockVoxel E [0 ; 19999]  -> ODFBlockVoxelIndex
voxel      E [0 ; 999999] -> Voxelindex

==> NewODFSampledImageArray[ voxel*642 + dir ] = ODFBlocks[ N ][ BlockVoxel*321 + dir ] ;
*/
  unsigned int BlockIndex = 0;
  unsigned int BlockVoxelIndex = 0;
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

    if( fa[ VoxelIndex ] == 0.0 && MaskArray[ VoxelIndex ] == 0.0 ) // if outside mask -> all ODFs = 0 // Outside of mask == both fa0 and index0 matrices show 0
    {
      for( unsigned int dir=0; dir < NBDIRS ; ++dir ) // for all directions
      {
        NewODFSampledImageArray[ VoxelIndex*NBDIRS + dir ] = 0.0;
      }
    } // if outside mask
    else // if inside mask
    {
      for( unsigned int dir=0; dir < NBDIRS ; ++dir ) // for all directions
      {
        NewODFSampledImageArray[ VoxelIndex*NBDIRS + dir ] = ODFBlocks[ BlockIndex ][ BlockVoxelIndex*NBDIRS + dir ];
      } // for all directions

      ++BlockVoxelIndex;
      if( BlockVoxelIndex >= BlockSizes[ BlockIndex ] ) // end of block: 20000 => last block will not go here but ok because no need to increment BlockIndex
      {
        BlockVoxelIndex = 0;
        ++BlockIndex;
      } // end of block
    } // if inside mask
  } // for all voxels in the image
  std::cout<< "100%" <<std::endl;

  // Write out sampled ODF image
  if( outputSampledODF != "" )
  {
    WriteOutITKODFImage( NewODFSampledImage, outputSampledODF );
  }

  return NewODFSampledImage; // NewODFSampledImage has 321 components
}

/// Master Function
ODFImageType::Pointer writeITK( std::string ODFfib, std::string outputCoeffsODFITK, std::string outputSampledODF ) // returns the assembled sampled ODF ITK image
{
  std::cout<<"Status : Reading fib file: " << ODFfib << std::endl;

  std::string outputFolder = itksys::SystemTools::GetRealPath( itksys::SystemTools::GetFilenamePath(outputCoeffsODFITK).c_str() );

  MatFile mat_reader;
  mat_reader.load_from_file( ODFfib.c_str() );
  if( mat_reader.get_matrix_count() == 0)
  {
    std::cout<<"Error  : Fib file empty" << std::endl;
    return NULL;
  }

  //// Get Data
  unsigned int row,col;

  // dim & spacing
  short* dim = 0;
  mat_reader.get_matrix("dimension",row,col,dim);
  std::cout<< row << " " << col << " | " << dim[0] <<" " <<  dim[1] << " " << dim[2] << " dimension"  <<std::endl;
  itk::ImageBase< 3 >::SizeType size;
  size[0] = dim[0];
  size[1] = dim[1];
  size[2] = dim[2];

  double* vs = 0;
  mat_reader.get_matrix("voxel_size",row,col,vs);
  std::cout<< row << " " << col << " | " << vs[0] <<" " <<  vs[1] << " " << vs[2] << " voxel_size" <<std::endl;
  itk::ImageBase< 3 >::SpacingType spacing;
  spacing[0] = vs[0];
  spacing[1] = vs[1];
  spacing[2] = vs[2];

  // odf vertices
  float* odfVertices = 0;
  mat_reader.get_matrix("odf_vertices",row,col,odfVertices);
  std::cout<< row << " " << col << " odf_vertices" <<std::endl;
  GetVertices< float >(odfVertices, col, (outputFolder + "/odf_vertices.txt").c_str());

  // odf faces
  short* odfFaces = 0;
  mat_reader.get_matrix("odf_faces",row,col,odfFaces);
  std::cout<< row << " " << col << " odf_faces" <<std::endl;
  GetVertices< short >(odfFaces, col, (outputFolder + "/odf_faces.txt").c_str());

  // Mask (index0)
  short* MaskArray;
  mat_reader.get_matrix("index0",row,col,MaskArray);
  std::cout<< row << " " << col << " index0" <<std::endl;
  WriteITKScalarImage< short >( MaskArray, outputFolder + "/index0.nrrd", size, spacing );

  // fa
  double* fa;
  mat_reader.get_matrix("fa0",row,col,fa);
  std::cout<< row << " " << col << " fa0" <<std::endl;
  WriteITKScalarImage< double >( fa, outputFolder + "/fa0.nrrd", size, spacing );

  // odf
  std::vector< ODFType* > ODFBlocks;
  std::vector< unsigned int > BlockSizes;
  for (unsigned int index = 0; index < mat_reader.get_matrix_count() ; ++index)
  {
    std::string MatrixName ( mat_reader.get_matrix_name(index) ) ;
//    std::cout<< MatrixName <<std::endl;
    if( MatrixName.find( (std::string)"odf" ) != std::string::npos  &&  MatrixName != "odf_vertices"  &&  MatrixName != "odf_faces" )
    {
      ODFType* odf;
      mat_reader.get_matrix( MatrixName.c_str(), row, col, odf );
      ODFBlocks.push_back( odf );
      BlockSizes.push_back( col );

      std::cout<< row << " " << col << " " << MatrixName <<std::endl;
    }
  }
  ODFImageType::Pointer NewODFSampledImage = AssembleODFBlocks( ODFBlocks, BlockSizes, outputSampledODF, size, spacing, fa, MaskArray ); // FA needed for masking
  if( ! NewODFSampledImage )
  {
    return NULL;
  }

  // Reconstruct sph. harm. from ODf samples
  ComputeSpharmCoeffs( NewODFSampledImage, outputCoeffsODFITK, size, spacing ); // NewODFSampledImage has 321 components

  return NewODFSampledImage;

} // bool writeITK()

  /////////////////////////////////////////
 //                MAIN                 //
/////////////////////////////////////////

int main (int argc, char* argv[])
{
  PARSE_ARGS;

/* Check input args */
  if( ODFfib == "" &&  ODFitk == "" )
  {
    std::cout<<"Error  : Please give an input ODF image."<<std::endl;
    return 1;
  }

  if( ODFfib != "" &&  ODFitk != "" )
  {
    std::cout<<"Warning: As both ITK and fib files have been given, the ITK image will be used."<<std::endl;
    ODFfib = ""; // so the ITK image is used
  }

  if( ODFitk != "" &&  fa == "" ) // TODO: FA can be any scalar image to display with the ODF in DSI studio
  {
    std::cout<<"Error  : A FA image is required to convert ITK image to fib file."<<std::endl;
    return 1;
  }

  if( ODFitk != "" &&  mask == "" ) // TODO: If no mask given, create a mask image of the size of the ODF filled with all ones
  {
    std::cout<<"Error  : A mask is required to convert ITK image to fib file."<<std::endl;
    return 1;
  }

  if( outputODF == "" )
  {
    if( ODFfib != "" )
    {
      outputODF = "./ODFcoeffs.nrrd";
    }
    else
    {
      outputODF = "./ODF.fib.gz";
    }
    std::cout<<"Warning: No output image given, the output image will be: "<< outputODF <<std::endl;
  }

/* Convert image */
  if( ODFfib != "" ) // fib -> ITK
  {
    ODFImageType::Pointer NewODFSampledImage = writeITK( ODFfib, outputODF, outputSampledODF );

    if( ! NewODFSampledImage )
    {
      return 1;
    }

  }
  else // ITK -> fib
  {
    ODFImageType::Pointer ODFSampledImage = loadITKODF( ODFitk, outputSampledODF );

    if( ! ODFSampledImage )
    {
      return 1;
    }

    if( ! writeFib( ODFSampledImage, outputODF, mask, fa ) )
    {
      return 1;
    }
  }

  return 0;
} // int main

/*

// Mask
src/libs/dsi/dsi_interface_imp.cpp:322 : mask reconstruction

src/libs/dsi/odf_process.hpp:74 : creation of odf_index_map[]


// WRITE .fib file

src/libs/dsi/odf_process.hpp
    -> 120: virtual void end(Voxel& voxel,MatFile& mat_writer) ->  mat_writer.add_matrix()
    -> 60: const unsigned int odf_block_size = 20000;
    -> struct OutputODF
    -> 75: OutputODF::init() -> ODF contain only voxels in mask
    -> std::vector<std::vector<float> > odf_data;  => ODF data stored in float

// .fib file

src/mainwindow.cpp: void MainWindow::loadFib(QString filename)
    -> std::auto_ptr<ODFModel> new_handle(new ODFModel);  -> new_handle = ptr on ODFModel
    -> new_handle->load_from_file( filename )

src/libs/tracking/tracking_model.hpp: class ODFModel
     -> ODFModel::load_from_file() -> fib_data.load_from_file()
     -> public:  FibData fib_data;

src/libs/tracking/fib_data.hpp: class FibData
     -> FibData::load_from_file() -> mat_reader.load_from_file() -> fib.add_data(mat_reader)
     -> public:  MatFile mat_reader;

src/libs/mat_file.hpp: class MatFile
     ===> MatFile::load_from_file()
     ===> MatFile::write_to_file()

src/libs/tracking/fib_data.hpp: class FiberDirection
     ===> FiberDirection::add_data()
     ===> FiberDirection contains ODFData

src/libs/tracking/fib_data.hpp: struct ODFData
     ===> ODFData contains odfs


// ODF vertices

src/libs/dsi/tessellated_icosahedron.hpp: class tessellated_icosahedron


// Name of .fib file

src/reconstruction/reconstruction_window.cpp: void reconstruction_window::doReconstruction()
    -> const char* msg = (const char*)reconstruction();

src/libs/dsi/dsi_interface_imp.cpp: const char* reconstruction()


// .src file

src/mainwindow.cpp: void MainWindow::loadSrc(QStringList filenames)
    ->  reconstruction_window* new_mdi = new reconstruction_window(filenames,this);

src/reconstruction/reconstruction_window.cpp: reconstruction_window::reconstruction_window()
    -> load_src(0);
    -> bool reconstruction_window::load_src(int index)
    -> handle->load_from_file(filenames[index].toLocal8Bit().begin())

src/reconstruction/reconstruction_window.h: class reconstruction_window:
    -> private:  std::auto_ptr<ImageModel> handle;

src/libs/dsi/image_model.hpp: struct ImageModel
    -> ImageModel.load_from_file() -> mat_reader->load_from_file()
    -> public:   std::auto_ptr<MatFile> mat_reader;

src/libs/mat_file.hpp: class MatFile
     ===> MatFile::load_from_file()
     ===> MatFile::write_to_file()

*/
