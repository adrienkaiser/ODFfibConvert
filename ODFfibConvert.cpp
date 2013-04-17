/*
./ODFfibConvert --ODFfib "/home/akaiser/Networking/NFG/nfg_1.1.1/generated_collections/130325112527/DWI/dwi-00.src.gz.odf8.f5rec.qbi.sh8.0.006.fib.gz"
*/

/* std classes */
#include <vector>
#include <fstream>
#include <iostream>

#include "zlib.h"
#include "ODFfibConvertCLP.h" //generated when ccmake

/* ITK classes */
#include <itkImage.h>
#include <itkImageBase.h> // for SizeType and SpacingType without having to declare an ImageType
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itksys/SystemTools.hxx>

#define DOUBLE  0
#define FLOAT   10
#define UINT    20
#define SHORT   30
#define USHORT  40
#define UCHAR   50

/*
char              1byte
short int (short) 2bytes
int               4bytes
long int (long)   4bytes
bool              1byte
float             4bytes
double            8bytes
long double       8bytes
*/
const unsigned int element_size_array[10] = {8,4,4,2,2,1,0,0,0,0};

  /////////////////////////////////////////
 //           MAT MATRIX CLASS          //
/////////////////////////////////////////

class MatMatrix
{
public:
  std::string name;
  unsigned int type;
  unsigned int rows;
  unsigned int cols;
  unsigned int count;
  unsigned int namelen;
  std::vector<char> data_buf;
  char* data_ptr;

  bool readMat(void* in)
  { 
    unsigned int imagf = 0; // content useless but variable necessary

    if (gzread(in,(char*)&type,4) == -1 || type > 100 || type % 10 > 1) return false;
    if (type % 10) type = 0; // text
    if (gzread(in,(char*)&rows,4) == -1) return false;
    if (gzread(in,(char*)&cols,4) == -1) return false;
    if (gzread(in,(char*)&imagf,4) == -1) return false;
    if (gzread(in,(char*)&namelen,4) == -1) return false;
    std::vector<char> buffer(namelen+1);
    if (gzread(in,(char*)&*buffer.begin(),namelen) == -1) return false;

    count = rows*cols;
    name = &*buffer.begin();

    unsigned int total_size = count*element_size_array[(type%100)/10];
    try
    {
      std::vector<char> allocator(total_size);
      allocator.swap(data_buf); // swap ??
    }
    catch (...)
    {
      return false;
    }
    data_ptr = &*data_buf.begin(); // data_buf only used to allocate memory

    return gzread(in,(char*)data_ptr,total_size) != -1;
  }

}; // class MatMatrix

  /////////////////////////////////////////
 //                 FIB                 //
/////////////////////////////////////////

std::vector<MatMatrix*> loadFib(std::string filename)
{
  std::vector<MatMatrix*> ImageMatrices;
  std::vector<MatMatrix*> EmptyVector; // EXIT_FAILURE

  std::cout<<"Status : Loading fib file: "<<filename<<std::endl;

  void* in = gzopen(filename.c_str(), "rb");
  if (!in)
  {
     std::cout<<"Error  : gzopen failed"<<std::endl;
     return EmptyVector;
  }

  for(unsigned int index = 0 ; !gzeof(in) ; ++index)
  {
    MatMatrix * matrix= new MatMatrix;

    if(!matrix->readMat(in))
    {
      std::cout<<"Error  : read matrix failed"<<std::endl;
      return EmptyVector;
    }

    ImageMatrices.push_back(matrix);

  } // for()

  gzclose(in);

  std::cout<<"Status : Done loading fib file"<<std::endl;

  return ImageMatrices;

} // std::vector<MatMatrix*> loadFib()

  /////////////////////////////////////////
 //                 ITK                 //
/////////////////////////////////////////

template<class T > void GetSize( char* data_ptr , itk::ImageBase< 3 >::SizeType &size )
{
  for( int i = 0 ; i < 3 ; i++ )
  {
    size[ i ] = static_cast< itk::ImageBase< 3 >::SizeValueType > ( *(T *)( data_ptr + sizeof(T) * i ) );  // char= 1 byte => get N bytes of the char array for each T
    //                                            SizeValueType = unsigned long
  }
}

template<class T > void GetSpacing( char* data_ptr , itk::ImageBase< 3 >::SpacingType &spacing )
{
  for( int i = 0 ; i < 3 ; i++ )
  {
    spacing[ i ] = static_cast< itk::ImageBase< 3 >::SpacingValueType > ( *(T *)( data_ptr + sizeof(T) * i ) );  // char= 1 byte => get N bytes of the char array for each T
    //                                               SpacingValueType = double
  }
}

template<class T > void WriteITKImage( char* data_ptr, std::string filename, typename itk::Image <T, 3>::SizeType size, typename itk::Image <T, 3>::SpacingType spacing )
{
  std::cout<<"Status : Writing out ITK image: " << filename << "...";

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

  // direction
//  ImageType::DirectionType direction;
//  NewImage->SetDirection ( direction ) ;

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

  // Cast data and fill buffer
  for( int i = 0 ; i < size[0]*size[1]*size[2] ; i++ )
  {
    NewImageArray[ i ] = *(T *)( data_ptr + sizeof(T) * i ) ;  // char= 1 byte => get N bytes of the char array for each T
  }

  // Write out image
  typedef itk::ImageFileWriter < ImageType > WriterType ;
  typename WriterType::Pointer Writer = WriterType::New() ;
  Writer->SetFileName ( filename ); 
  Writer->SetInput( NewImage );
  Writer->SetUseCompression(true);

  try
  {
    Writer->Update() ;
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Problem writing the image: " << filename << std::endl ;
    std::cerr << excp << std::endl ;
  }

  std::cout<<"DONE"<<std::endl;

  return;
} // void WriteITKImage()


template<class T > void GetVertices( char* data_ptr, int nbVertices, std::string filename ) // matrix dim = (3, nbVertices)
{
  std::cout<<"Status : Writing file: " << filename <<std::endl;

  std::vector< std::vector< T > > Vertices;
  Vertices.resize(3);

  // Get vertices
  for( int i = 0 ; i < 3 ; i++ )
  {
    Vertices[i].resize(nbVertices);
    for( int j = 0 ; j < nbVertices ; j++ )
    {
      Vertices[i][j] = *(T *)( data_ptr + sizeof(T) * i*j ) ;  // char= 1 byte => get N bytes of the char array for each T
    }
  }

  // Write out in a file
  std::ofstream VerticesFileStream (filename.c_str() , std::ios::out | std::ios::trunc);  // opening in writing with erasing the open file
  if(! VerticesFileStream)
  {
    std::cout<<"Error  : Creating file: "<< filename <<std::endl;
    return;
  }

  for( int j = 0 ; j < nbVertices ; j++ ) // 1 vertex by line
  {
    VerticesFileStream << Vertices[0][j] << " " << Vertices[1][j] << " " << Vertices[2][j] << std::endl ;
  }

  VerticesFileStream.close();

  return;
}

bool writeITK( std::vector<MatMatrix*> ImageMatrices, std::string outputODF)
{
  std::cout<<"Status : Writing ITK image: "<<outputODF<<std::endl;

  std::string outputFolder = itksys::SystemTools::GetRealPath( itksys::SystemTools::GetFilenamePath(outputODF).c_str() );

  itk::ImageBase< 3 >::SizeType size;
  itk::ImageBase< 3 >::SpacingType spacing;
  float z0;

  std::vector< std::string > ImagesToWriteOut; // = {"gfa", "iso", "fa0", "fa1", "fa2", "fa3", "fa4", "nqa0", "nqa1", "nqa2", "nqa3", "nqa4", "index0", "index1", "index2", "index3", "index4"};
  ImagesToWriteOut.push_back("gfa");
  ImagesToWriteOut.push_back("iso");
  ImagesToWriteOut.push_back("fa0");
  ImagesToWriteOut.push_back("fa1");
  ImagesToWriteOut.push_back("fa2");
  ImagesToWriteOut.push_back("fa3");
  ImagesToWriteOut.push_back("fa4");
  ImagesToWriteOut.push_back("nqa0");
  ImagesToWriteOut.push_back("nqa1");
  ImagesToWriteOut.push_back("nqa2");
  ImagesToWriteOut.push_back("nqa3");
  ImagesToWriteOut.push_back("nqa4");
  ImagesToWriteOut.push_back("index0");
  ImagesToWriteOut.push_back("index1");
  ImagesToWriteOut.push_back("index2");
  ImagesToWriteOut.push_back("index3");
  ImagesToWriteOut.push_back("index4");

  // for all matrices
  for (unsigned int index = 0; index < ImageMatrices.size() ;++index)
  {
    // Dimension
    if( ImageMatrices[index]->name == "dimension" )
    {
      switch( ImageMatrices[index]->type )
      {
        case UINT:
          GetSize< unsigned int >(ImageMatrices[index]->data_ptr, size);
          break ;
        case SHORT:
          GetSize< short >(ImageMatrices[index]->data_ptr, size);
          break ;
        case USHORT:
          GetSize< unsigned short >(ImageMatrices[index]->data_ptr, size);
          break ;
        case UCHAR:
          GetSize< unsigned char >(ImageMatrices[index]->data_ptr, size);
          break ;
        default:
          return false;
      } // switch()

      std::cout<<"Dimension = ["<< size[0] << ", "<< size[1]<< " , "<< size[2]<<"]"<<std::endl; 
    } // if dimension

    // Spacing
    else if( ImageMatrices[index]->name == "voxel_size" )
    {
      switch( ImageMatrices[index]->type )
      {
        case DOUBLE:
          GetSpacing< double >(ImageMatrices[index]->data_ptr, spacing);
          break ;
        case FLOAT:
          GetSpacing< float >(ImageMatrices[index]->data_ptr, spacing);
          break ;
        default:
          return false;
      }
      std::cout<<"VoxelSize = ["<< spacing[0] << ", "<< spacing[1]<< " , "<< spacing[2]<<"]"<<std::endl;
    } // if spacing

    // z0
    else if( ImageMatrices[index]->name == "z0" )
    {
      z0 = *(float *)ImageMatrices[index]->data_ptr ;
      std::cout<<"z0        = "<<z0<<std::endl;
    }

    // images to write out
//    else if( ImageMatrices[index]->name == "gfa" ) // !!! any image needs to be after size & spacing
    // if vector contains name of current matrix -> write out image
    else if( std::find(ImagesToWriteOut.begin(), ImagesToWriteOut.end(), (std::string)ImageMatrices[index]->name) != ImagesToWriteOut.end() )
    {
      switch( ImageMatrices[index]->type )
      {
        case DOUBLE:
          WriteITKImage< double >(ImageMatrices[index]->data_ptr, outputFolder + "/" + ImageMatrices[index]->name + ".nrrd", size, spacing);
          break ;
        case FLOAT:
          WriteITKImage< float >(ImageMatrices[index]->data_ptr, outputFolder + "/" + ImageMatrices[index]->name + ".nrrd", size, spacing);
          break ;
        case SHORT:
          WriteITKImage< short >(ImageMatrices[index]->data_ptr, outputFolder + "/" + ImageMatrices[index]->name + ".nrrd", size, spacing);
          break ;
        case UINT:
          WriteITKImage< unsigned int >(ImageMatrices[index]->data_ptr, outputFolder + "/" + ImageMatrices[index]->name + ".nrrd", size, spacing);
          break ;
        default:
          return false;
      }
    } // if image to write out

    // odf_vertices and odf_faces
    else if( ImageMatrices[index]->name == "odf_vertices" || ImageMatrices[index]->name == "odf_faces")
    {
      switch( ImageMatrices[index]->type )
      {
        case DOUBLE:
          GetVertices< double >(ImageMatrices[index]->data_ptr, ImageMatrices[index]->cols, outputFolder + "/" + ImageMatrices[index]->name + ".txt");
          break ;
        case FLOAT:
          GetVertices< float >(ImageMatrices[index]->data_ptr, ImageMatrices[index]->cols, outputFolder + "/" + ImageMatrices[index]->name + ".txt");
          break ;
        case SHORT:
          GetVertices< short >(ImageMatrices[index]->data_ptr, ImageMatrices[index]->cols, outputFolder + "/" + ImageMatrices[index]->name + ".txt");
          break ;
        case UINT:
          GetVertices< unsigned int >(ImageMatrices[index]->data_ptr, ImageMatrices[index]->cols, outputFolder + "/" + ImageMatrices[index]->name + ".txt");
          break ;
        default:
          return false;
      }
    } // if odf_vertices or odf_faces

  } // for all matrices


  return true;

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
  }

  if( ODFfib != "" &&  ODFitk != "" )
  {
    std::cout<<"Warning: As both ITK and fib files have been given, the ITK image will be used."<<std::endl;
    ODFfib = ""; // so the ITK image is used
  }

  if( outputODF == "" )
  {
    if( ODFfib != "" )
    {
      outputODF = "./ODF.nrrd";
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
    writeITK( loadFib(ODFfib), outputODF );
  }
  else // ITK -> fib
  {
    //writeFib( loadITK(ODFitk), outputODF );
  }

  return 0;

} // main

/* Display matrices content

    std::cout<<"======== MAT ========================"<<std::endl;

    std::cout<<"type\trows\tcols\tnamelen\tcount\tname"<<std::endl;
    std::cout<<ImageMatrices[index]->type<<"\t";
    std::cout<<ImageMatrices[index]->rows<<"\t";
    std::cout<<ImageMatrices[index]->cols<<"\t";
    std::cout<<ImageMatrices[index]->namelen<<"\t";
    std::cout<<ImageMatrices[index]->count<<"\t";
    std::cout<<ImageMatrices[index]->name<<"\t"<<std::endl;
    if( ImageMatrices[index]->count<100 ) 
    {
      std::cout<<"===================================="<<std::endl;
      std::cout<<"data_ptr=|"<< ImageMatrices[index]->data_ptr <<"|"<<std::endl;
      std::cout<<"===================================="<<std::endl;
    }

*/

/*

// WRITE .fib file

src/libs/dsi/odf_process.hpp:120 -> virtual void end(Voxel& voxel,MatFile& mat_writer) ->  mat_writer.add_matrix()


// .fib file

src/mainwindow.cpp: void MainWindow::loadFib(QString filename)
    -> std::auto_ptr<ODFModel> new_handle(new ODFModel);  -> new_handle = ptr on ODFModel
    -> new_handle->load_from_file( filename )

src/libs/tracking/tracking_model.hpp: class ODFModel
     -> ODFModel::load_from_file() -> fib_data.load_from_file()
     -> public:  FibData fib_data;

src/libs/tracking/fib_data.hpp: class FibData
     -> FibData::load_from_file() -> mat_reader.load_from_file()
     -> public:  MatFile mat_reader;

src/libs/mat_file.hpp: class MatFile
     ===> MatFile::load_from_file()
     ===> MatFile::write_to_file()

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
