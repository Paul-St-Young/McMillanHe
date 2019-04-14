#pragma once

#include <iostream>
#include "hdf5.h"
#include <Eigen/Dense>
#include <vector>

namespace ezh5{
  ///////////////////////////
  class ID{
    public:
    hid_t id;
    ID(): id(-1) {}
    ID(hid_t id_in) : id(id_in) {}
    ~ID(){}
  };
  ///////////////////////////
  ///HDF5 nodes
  class Node :public ID
  {
  public:
    hid_t pid,DataSpaceID;
    Eigen::Matrix<hsize_t,Eigen::Dynamic,1> dims;
    std::string path;
    bool verbose_;

    Node(){verbose_ = false;}
    Node(hid_t pid_in, const std::string& path_in, bool verbose=false): ID(-1), pid(pid_in), DataSpaceID(-1),path(path_in){
      verbose_ = verbose;
      if(verbose_) std::cout<<"creating "<<path<<std::endl;
    }
    ////////
    Node& operator()(const std::string& path){
      return *this;
    }
    void flush(){ H5Fflush(this->id,H5F_SCOPE_GLOBAL);}
    ////////
    Node operator[] (const std::string& path_more);
    ////////
    template<typename T> Node& operator =  (T val);
    template<typename T> Node& operator >> (T& val);
    template<typename T> Node& operator << (T val);
    ////////
    template<typename T> Node& operator =  (std::vector<T>& vec);
    template<typename T> Node& operator >> (std::vector<T>& vec);
    template<typename T> Node& operator << (std::vector<T>& vec);
    ////////
    template<typename T> Node& operator =  (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor>& m);
    template<typename T> Node& operator >> (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor>& m);

    template<typename T> Node& operator =  (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor>& m);
    template<typename T> Node& operator >> (Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor>& m);

    template<typename T> Node& operator =  (Eigen::Matrix<T, Eigen::Dynamic, 1>& m);//vector
    template<typename T> Node& operator >> (Eigen::Matrix<T, Eigen::Dynamic, 1>& m);//vector
    template<typename T> Node& operator << (Eigen::Matrix<T, Eigen::Dynamic, 1>& m);//vector
    ////////
    template<typename T> Node& append(T val);
    template<typename T> Node& append(std::vector<T>& vec);
    template<typename T> Node& append(Eigen::Matrix<T, Eigen::Dynamic, 1>& m);
    ////////
    ~Node(){
      if(verbose_) std::cout<<"closing "<<path<<std::endl;
      if (this->id>0 && this->DataSpaceID>0){
        H5Dclose(this->id);
        this->id=-1;
        H5Sclose(this->DataSpaceID);
        this->DataSpaceID=-1;
      }
      else if(this->id>0)
      {
        H5Gclose(this->id);
        this->id = -1;
      }
    }
  };
  ///////////////////////////
  /** @brief an HDF5 File
   * @details a simple interface for reading/writing to an HDF5 file
   * See the following example for usage
   *    {
        //put inside brackets so the deconstructor is called
        ezh5::File f ("universe.h5", H5F_ACC_TRUNC);  // create a file object called f
        f["i"]=1;
        f["d"]=5.2;
        f["c"]=1.0+5.1*CONE;
        f["ch"]='t';
        std::vector<double> v(10);
        double x=0.0;
        for(auto &vi : v){
          vi=x;
          x+=2.5;
        }
        f["v"]=v;
        Eigen::VectorXd ev(10);
        ev.setRandom();
        f["ev"]=ev;
        Eigen::VectorXcd evc(10);
        evc.setRandom();
        f["evc"]=evc;
        Eigen::MatrixXd mxd(5,5);
        mxd.setRandom();
        f["mxd"]=mxd;
        
        Eigen::MatrixXcd mxcd(5,5);
        mxcd.setRandom();
        f["mxcd"]=mxcd;
        
        }

        {
        ezh5::File g ("universe.h5", H5F_ACC_RDONLY);  // create a file object called f
        int i;
        g["i"] >> i;
        double d;
        g["d"] >> d;
        std::complex<double> c;
        g["c"] >> c;
        char name;
        g["ch"] >> name;
        std::vector<double> v;
        g["v"] >> v;
        Eigen::VectorXd ev;
        g["ev"] >> ev;
        Eigen::VectorXcd evc;
        g["evc"] >> evc;
        Eigen::MatrixXd mxd;
        g["mxd"] >> mxd;
        Eigen::MatrixXcd mxcd;
        g["mxcd"] >> mxcd;

        std::cout<<"i:"<<i<<std::endl;
        std::cout<<"d:"<<d<<std::endl;
        std::cout<<"c:"<<c<<std::endl;
        std::cout<<"ch:"<<name<<std::endl;
        std::cout<<"v:";
        for(auto x: v)
          std::cout<<x<<" ";
        std::cout<<std::endl;
        std::cout<<"ev:"<<ev.transpose()<<std::endl;  
        std::cout<<"evc:"<<evc.transpose()<<std::endl;
        std::cout<<"mxd:"<<mxd<<std::endl;
        std::cout<<"mxcd:"<<mxcd<<std::endl;
        g["mxd"] >> mxcd;
        std::cout<<"mxcd:"<<mxcd<<std::endl;
        }

   */
  class File : public Node
  {
  public:
    File(bool verbose = false){verbose_=verbose;this->id=-1;};
    File(const std::string& path, unsigned flags, bool verbose = false);
    void open(const std::string& path, unsigned flags, bool verbose = false);

    ~File()
    {
      if(this->id != -1) H5Fclose(this->id);
      this->id = -1;
    }
  };
  ///////////////////////////
}
