
#include <vector>
#include <fstream>
#include <iostream>
#include "zlib.h"

/*
0 : double
10: float
20: unsigned int
30: short
40: unsigned short
50: unsigned char
*/

const unsigned int element_size_array[10] = {8,4,4,2,2,1,0,0,0,0};

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


bool load_from_file(const char* file_name)
{
  std::string filename = file_name;

  void* in = gzopen(file_name, "rb");
  if (!in)
  {
     std::cout<<"gzopen failed"<<std::endl;
     return false;
  }

  for (unsigned int index = 0;!gzeof(in);++index)
  {
    std::cout<<"======== MAT ========================"<<std::endl;

    MatMatrix * matrix= new MatMatrix;

    if (!matrix->readMat(in))
    {
      std::cout<<"read matrix failed"<<std::endl;
      return false;
    }

    std::cout<<"type\trows\tcols\tnamelen\tcount\tname"<<std::endl;
    std::cout<<matrix->type<<"\t"<<matrix->rows<<"\t"<<matrix->cols<<"\t"<<matrix->namelen<<"\t"<<matrix->count<<"\t"<<matrix->name<<"\t"<<std::endl;
    if( matrix->count<100 ) 
    {
      std::cout<<"===================================="<<std::endl;
      std::cout<<"data_ptr=|"<< matrix->data_ptr <<"|"<<std::endl;
      std::cout<<"===================================="<<std::endl;
    }

  } // for()

  gzclose(in);

  return true;

} // bool load_from_file()


int main (int argc, char* argv[])
{
  //  load_from_file("/home/akaiser/Downloads/CMU_30_20130213build.gz.mean.fib.gz");
    load_from_file("/home/akaiser/Networking/NFG/nfg_1.1.1/generated_collections/130325112527/DWI/dwi-00.src.gz.odf8.f5rec.qbi.sh8.0.006.fib.gz");

    return 0;
}

/*

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
