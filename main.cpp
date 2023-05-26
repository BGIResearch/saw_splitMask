#include <iostream>
#include <string>
#include <fstream>
#include <thread>
#include <mutex>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <set>
#include <map>
#include <vector>
#include <stdlib.h>
#include <hdf5.h>
using namespace ::std;
#define DATASETNAME "bpMatrix_"
#define RANK 3
#define MAJOR_VERSION 2
#define MINOR_VERSION 0.1
enum errCode : unsigned int
{
  // 01-19
  err_sw_exception = 6,
  err_param_invalid=1,
  err_fileOpen_failed,
  err_fileParse_failed,
  err_fileOtherIO_failed,
  err_otherAPI_failed
};
void sawErrCode(unsigned int errCode, string message)
{
  string errCodeStr = to_string(errCode);
  int s=errCodeStr.size();
  if (errCodeStr.size() < 3)
  {
    for (string::size_type i = 0; i < 3 - s; i++)
    {
      errCodeStr = "0" + errCodeStr;
    }
  }
  time_t timep;
  time(&timep);
  struct tm *p = localtime(&timep);
  string outStr = "[" + to_string(p->tm_year + 1900) + "-" + to_string(p->tm_mon + 1) + "-" + to_string(p->tm_mday) + " " + to_string(p->tm_hour) + ":" + to_string(p->tm_min) + ":" + to_string(p->tm_sec) + "]\t" + "SAW-A00" + errCodeStr + ": " + message;
  ofstream errLog("errcode.log", ios::app);
  errLog << outStr << endl;
  errLog.close();
}
void printVersionAndAuthor(){
  cout<<"Program: splitBcBin\n";
  cout<<"Version: "<<MAJOR_VERSION<<"."<<MINOR_VERSION<<endl;
//   cout<<"Contact: GongChun<gongchun@genomics.cn>"<<endl;
}
const uint8_t RC_BASE[4] = {2, 3, 0, 1};
const uint8_t RC_BASE2[4] = {3, 2, 1, 0};
const char ATCG_BASES[4] = {'A', 'C', 'T', 'G'};
const char ATCG_BASES2[4] = {'A', 'C', 'G', 'T'};
uint64_t old2new[4] = {0, 1, 3, 2};
// A C G T = 0 1 2 3
const short base2int[24] = {99, 0, 99, 1, 99, 99, 99, 2, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 3, 99, 99, 99};
uint64_t fileSize = 0;
int threadsNum = 16;
int splitNum = 64;
int fenmu = (1 << 16) / splitNum;
string barcodePositionMapFile;
string outDir;
int barcodeIntByteSize = 8;
int positionByteSize = 8;
// int barcodeLenInBin=25;
int outputBarcodeLen = 24;
int barcodeClassifyNum = 8;
string splitBcPos = "2_25";
void *subThread(int threadIdx);
void *subThreadForH5(int threadIdx);
mutex loglock;
int *fds;
// hid_t datasetID;
typedef struct Position1
{
  //	friend class boost::serialization::access;
  uint32_t x;
  uint32_t y;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &x;
    ar &y;
  }
} Position1;
uint64_t seqEncode(const char *sequence, const int &seqStart, const int &seqLen, bool isRC = false)
{

  uint64_t n = 0;
  uint64_t k = 0;
  for (int i = seqStart; i < seqLen; i++)
  {
    n = (sequence[i] & 6) >> 1; // 6:  ob00000110
    if (isRC)
    {
      n = RC_BASE[n];
      k |= (n << ((seqLen - i - 1) * 2));
    }
    else
    {
      k |= (n << (i * 2));
    }
  }
  return k;
}
string seqDecode(const uint64_t &seqInt, const int &seqLen)
{
  uint8_t tint;
  string seqs = "";
  for (int i = 0; i < seqLen; i++)
  {
    tint = (seqInt >> (i * 2)) & 3;
    seqs.push_back(ATCG_BASES[tint]);
  }
  return seqs;
}
uint64_t seqEncode2(const char *sequence, const int &seqStart, const int &seqLen, bool isRC = false)
{

  uint64_t n = 0;
  uint64_t k = 0;
  for (int i = seqLen - 1; i >= seqStart; i--)
  {
    //		n = (sequence[i] & 6) >> 1; //6:  ob00000110
    n = base2int[sequence[i] - 64];
    if (isRC)
    {
      n = RC_BASE2[n];
      k |= (n << ((seqLen - i - 1) * 2));
    }
    else
    {
      k |= (n << (i * 2));
    }
  }
  return k;
}
string seqDecode2(const uint64_t &seqInt, const int &seqLen)
{
  uint8_t tint;
  string seqs = "";
  for (int i = 0; i < seqLen; i++)
  {
    tint = (seqInt >> ((seqLen - i - 1) * 2)) & 3;
    seqs.push_back(ATCG_BASES2[tint]);
  }
  return seqs;
}
int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    printVersionAndAuthor();
    cerr << argv[0] << "\t<barcodefile.bin/.h5> <outputDir> <threads number|optional>  <splits number|optional> <barcode position: 1_24 or 2_25 |optional>" << endl;
    sawErrCode(err_param_invalid, "parameters error");
    exit(1);
  }
  barcodePositionMapFile = argv[1];
  outDir = argv[2];
  if (argc > 3)
  {
    threadsNum = atoi(argv[3]);
  }
  if (argc > 4)
  {
    splitNum = atoi(argv[4]);
    fenmu = (1 << 16) / splitNum;
  }
  //    if(argc>5){
  //        barcodeLenInBin=atoi(argv[5]);
  //    }
  if (argc > 5)
  {
    splitBcPos = argv[5];
  }
  thread t_array[threadsNum];
  fds = new int[splitNum];
  string relaname = barcodePositionMapFile;
  if (barcodePositionMapFile.find("/") != string::npos)
  {
    relaname = barcodePositionMapFile.substr(barcodePositionMapFile.find_last_of("/") + 1);
  }
  if (splitNum >= 100)
  {
    cerr << "Error, splitNum is too large" << endl;
    // cerr << "please contact the author(gongchun@genomics.cn) to update the program" << endl;
    string errMsg="splitNum is too large";
    sawErrCode(err_param_invalid,errMsg);
    exit(1);
  }
  for (int i = 0; i < splitNum; i++)
  {
    string pad0 = "";
    if ((i + 1) < 10)
    {
      pad0 = "0";
    }
    string outfilename = outDir + "/" + pad0 + to_string(i + 1) + "." + relaname;
    if (outfilename.rfind(".h5") == outfilename.size() - 3)
    {
      outfilename = outfilename.replace(outfilename.rfind(".h5"), 3, ".bin");
    }
    fds[i] = open(outfilename.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0755);
    if (fds[i] < 0)
    {
      cerr << "Error,cannot write to file," << outfilename << endl;
      string errMsg="cannot write to file,"+outfilename;
      sawErrCode(err_fileOtherIO_failed,errMsg);
      exit(1);
    }
  }
  if (barcodePositionMapFile.rfind(".bin") == barcodePositionMapFile.size() - 4)
  {
    ifstream mapReader(barcodePositionMapFile, ios::in | ios::binary);
    if (!mapReader.is_open())
    {
      throw invalid_argument("Could not open the file: " + barcodePositionMapFile);
    }
    mapReader.seekg(0, ios::end);
    fileSize = mapReader.tellg();
    if (fileSize % (barcodeIntByteSize + positionByteSize) != 0L)
    {
      cout << argv[1] << "\tfileSize: " << fileSize << "\t" << barcodeIntByteSize + positionByteSize << endl;
      cerr << "Error, .bin file seems not right" << endl;
      string errMsg=".bin file seems not right";
      sawErrCode(err_fileParse_failed,errMsg);
      exit(1);
    }
    mapReader.close();

    for (int i = 0; i < threadsNum; i++)
    {
      t_array[i] = thread(subThread, i);
    }
    for (int i = 0; i < threadsNum; i++)
    {
      t_array[i].join();
    }
    for (int i = 0; i < splitNum; i++)
    {
      close(fds[i]);
    }
    //		boost::archive::binary_iarchive ia(mapReader);
    //		ia >> bpmap;
    // while (!mapReader.eof()) {
    //	mapReader.read((char*)&barcodeInt, sizeof(barcodeInt));
    //	mapReader.read((char*)&position, sizeof(position));
    //	bpmap[barcodeInt] = position;
    //}
  }
  else if (barcodePositionMapFile.rfind(".h5") == barcodePositionMapFile.size() - 3)
  {
    ifstream mapReader(barcodePositionMapFile, ios::in | ios::binary);
    if (!mapReader.is_open())
    {
      throw invalid_argument("Could not open the file: " + barcodePositionMapFile);
    }
    mapReader.close();

    for (int i = 0; i < threadsNum; i++)
    {
      t_array[i] = thread(subThreadForH5, i);
    }
    for (int i = 0; i < threadsNum; i++)
    {
      t_array[i].join();
    }
  }
  else
  {
    cerr << "Error,only support .bin or .h5 file" << endl;
    string errMsg="only support .bin or .h5 file";
    sawErrCode(err_fileParse_failed,errMsg);
    exit(1);
  }
  cout << "split task done" << endl;
  return 0;
}
void *subThread(int threadIdx)
{
  loglock.lock();
  uint64_t totalBcCount = fileSize / 16;
  uint64_t numEachThread = totalBcCount / threadsNum;
  uint64_t startPos = numEachThread * 16L * (uint64_t)threadIdx;
  uint64_t endPos = numEachThread * 16 * (threadIdx + 1);

  //为了将那些对应多个位置的barcode删除，不再文件分块读取，而是每个线程处理同一个拆分小类的数据
  if (threadIdx == threadsNum - 1)
  {
    endPos = fileSize;
  }

  if ((endPos - startPos) % (barcodeIntByteSize + positionByteSize) != 0)
  {
    cerr << "Error,code error" << endl;
    sawErrCode(err_sw_exception,"code error");
    exit(1);
  }
  cout << threadIdx << "\t" << startPos << "\t" << endPos << endl;
  //    cout<<"totalNum:\t"<<totalBcCount<<endl;
  loglock.unlock();
  //    return nullptr;
  ifstream mapReader(barcodePositionMapFile, ios::in | ios::binary);
  mapReader.seekg(startPos);
  uint64_t barcodeInt = 0;
  uint64_t coor = 0;

  // 前面8个碱基对应的数/1024+1
  int tagNum = 0;
  uint64_t iter = 0;
  uint64_t newBarcodeInt = 0;
  char *record = new char[16];
  //    map<uint64_t,int> bc1to24num,bc2to25num;
  map<uint64_t, vector<uint64_t>> bc1to24, bc2to25;
  //    bool a24=false;
  //    bool a25=false;
  //    ofstream outBc24("s2.dup24.bc.txt");
  //    ofstream outBc25("s2.dup25.bc.txt");
  while (mapReader.tellg() < endPos)
  {
    mapReader.read((char *)&barcodeInt, sizeof(barcodeInt));
    mapReader.read((char *)&coor, sizeof(coor));
    if (splitBcPos == "1_24")
    {
      for (int i = 0; i < outputBarcodeLen; i++)
      {
        if (i < barcodeClassifyNum)
        {
          tagNum |= old2new[(barcodeInt >> (i * 2)) & 3] << ((barcodeClassifyNum - i - 1) * 2);
        }
        newBarcodeInt |= (old2new[(barcodeInt >> (i * 2)) & 3] & 3) << ((outputBarcodeLen - i - 1) * 2);
      }
      //		  bc1to24[newBarcodeInt].emplace_back(coor);
    }
    else if (splitBcPos == "2_25")
    {
      for (int i = 1; i < outputBarcodeLen + 1; i++)
      {
        if (i - 1 < barcodeClassifyNum)
        {
          tagNum |= old2new[(barcodeInt >> (i * 2)) & 3] << ((barcodeClassifyNum - i) * 2);
        }
        newBarcodeInt |= (old2new[(barcodeInt >> (i * 2)) & 3] & 3) << ((outputBarcodeLen - i) * 2);
      }
      //		  bc2to25[newBarcodeInt].emplace_back(coor);
    }
    else
    {
      cerr << "Error, splitBcPos error" << endl;
      string errMsg="splitBcPos error, expected 1_24 or 2_25";
      sawErrCode(err_param_invalid,errMsg);
      exit(1);
    }

    int outPrefix = tagNum / fenmu;
    tagNum = 0;
    memcpy(record, (char *)&newBarcodeInt, sizeof(newBarcodeInt));
    memcpy(record + sizeof(newBarcodeInt), (char *)&coor, sizeof(coor));
    write(fds[outPrefix], record, 16);
    iter++;
    newBarcodeInt = 0L;

    //        break;
    if (iter % 1000000 == 0)
    {
      cout << "thread " << threadIdx << " processed: " << iter << endl;
    }
  }
  delete[] record;
  //    cout<<"totalNum:\t"<<totalNum<<"\t"<<totalBcCount<<endl;
  //    cout<<"1_24 bp size after rmdup:\t"<<exists_bc1to24.size()<<endl;
  //    cout<<"2_25 bp size after rmdup:\t"<<exists_bc2to25.size()<<endl;
}
void *subThreadForH5(int threadIdx)
{
  herr_t status;
  hid_t fileID = H5Fopen(barcodePositionMapFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  string datasetName = DATASETNAME + to_string(1);
  hid_t datasetID = H5Dopen(fileID, datasetName.c_str(), H5P_DEFAULT);

  hid_t dspaceID = H5Dget_space(datasetID);
  // rank=3
  int rank = H5Sget_simple_extent_ndims(dspaceID);
  hsize_t dims[rank];
  status = H5Sget_simple_extent_dims(dspaceID, dims, NULL);
  uint64_t matrixLen = 1;
  for (int i = 0; i < rank; i++)
  {
    //	cout<<i<<"\t"<<dims[i]<<endl;
    matrixLen *= dims[i];
  }
  bool processFlag = false;
  char *record = new char[16];
  //  int packInThread=dims[0]/threadsNum;
  //  int endPos=threadIdx==threadsNum-1?dims[0]:(threadIdx+1)*packInThread;
  int blockSize = dims[0] / (threadsNum * 2);
  int iter = 0;
  for (uint32_t r = 0; r < dims[0]; r += blockSize)
  {
    if ((r / blockSize) % threadsNum != threadIdx)
    {
      continue;
    }
    uint64_t cpSize = r + blockSize > dims[0] ? (uint64_t)(dims[0] - r) * dims[1] * dims[2] : (uint64_t)dims[1] * dims[2] * blockSize;
    //    cout<<"thread "<<threadIdx<<" : "<<r<<"\t"<<r+cpSize/(dims[1]*dims[2])<<"\t"<<dims[0]<<endl;
    uint64_t *bpMatrix_buffer = new uint64_t[cpSize];
    //	cout<<r<<"\t"<<blockSize<<"\t"<<cpSize<<endl;
    for (int i = 0; i < cpSize; i++)
    {
      bpMatrix_buffer[i] = 0L;
    }
    //	cout<<bpMatrix_buffer[0]<<endl;
    hsize_t dimsm2[1];
    dimsm2[0] = cpSize;
    hid_t memspace2 = H5Screate_simple(1, dimsm2, NULL);

    hsize_t offset[3], count[3];
    offset[0] = r;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = r + blockSize > dims[0] ? dims[0] - r : blockSize;
    count[1] = dims[1];
    count[2] = 1;
    hsize_t offset_out2[1];
    hsize_t count_out2[1];
    offset_out2[0] = 0;
    count_out2[0] = cpSize;
    if (count[0] * count[1] != cpSize)
    {
      cerr << "different size\t" << r << "\t" << count[0] << "\t" << count[1] << "\t" << cpSize << endl;
      sawErrCode(err_sw_exception,"code error");
      exit(1);
    }
    status = H5Sselect_hyperslab(memspace2, H5S_SELECT_SET, offset_out2, NULL,
                                 count_out2, NULL);
    status = H5Sselect_hyperslab(dspaceID, H5S_SELECT_SET, offset, NULL,
                                 count, NULL);
    status = H5Dread(datasetID, H5T_NATIVE_ULONG, memspace2, dspaceID,
                     H5P_DEFAULT, bpMatrix_buffer);
    int segment = 0;
    int tagNum = 0;
    uint64_t newBarcodeInt = 0L;

    uint64_t coor = 0;
    for (uint32_t i = 0; i < cpSize / (dims[1] * dims[2]); i++)
    {
      for (uint32_t c = 0; c < dims[1]; c++)
      {
        // bpMatrix[r][c] = bpMatrix_buffer + r*dims[1]*dims[2] + c*dims[2];
        //	  Position1 position = {c, r};
        if (rank >= 3)
        {
          segment = dims[2];
          for (int s = 0; s < segment; s++)
          {
            uint64_t barcodeInt = bpMatrix_buffer[i * dims[1] * dims[2] + c * segment + s];
            //			cout << rowIdx*dims[1]*dims[2]+c*segment + s << "\t" << (int)barcodeInt << endl;
            //			cout<<(int32_t)barcodeInt<<endl;
            if (barcodeInt == 0)
            {
              continue;
            }
            if (splitBcPos == "1_24")
            {
              for (int j = 0; j < outputBarcodeLen; j++)
              {
                if (j < barcodeClassifyNum)
                {
                  tagNum |= old2new[(barcodeInt >> (j * 2)) & 3] << ((barcodeClassifyNum - j - 1) * 2);
                }
                newBarcodeInt |= (old2new[(barcodeInt >> (j * 2)) & 3] & 3) << ((outputBarcodeLen - j - 1) * 2);
              }
            }
            else if (splitBcPos == "2_25")
            {
              for (int j = 1; j < outputBarcodeLen + 1; j++)
              {
                if (j - 1 < barcodeClassifyNum)
                {
                  tagNum |= old2new[(barcodeInt >> (j * 2)) & 3] << ((barcodeClassifyNum - j) * 2);
                }
                newBarcodeInt |= (old2new[(barcodeInt >> (j * 2)) & 3] & 3) << ((outputBarcodeLen - j) * 2);
              }
            }
            else
            {
              cerr << "Error, splitBcPos error" << endl;
              string errMsg="splitBcPos error, expected 1_24 or 2_25";
              sawErrCode(err_param_invalid,errMsg);
              exit(1);
            }

            int outPrefix = tagNum / fenmu;
            tagNum = 0;
            memcpy(record, (char *)&newBarcodeInt, sizeof(newBarcodeInt));
            memcpy(record + sizeof(newBarcodeInt), (char *)&c, sizeof(c));
            uint32_t realRow = r + i;
            memcpy(record + sizeof(newBarcodeInt) + sizeof(c), (char *)&realRow, sizeof(realRow));
            write(fds[outPrefix], record, 16);
            //			write(smallFds[outPrefix], record, 16);
            iter++;
            newBarcodeInt = 0L;

            //        break;
            if (iter % 1000000 == 0)
            {
              cout << "thread " << threadIdx << " processed: " << iter << endl;
            }
            //		  bpMap[barcodeInt] = position;
          }
        }
        else
        {
          uint64_t barcodeInt = bpMatrix_buffer[c];
          if (barcodeInt == 0)
          {
            continue;
          }
          //		bpMap[barcodeInt] = position;
        }
      }
    }
    delete[] bpMatrix_buffer;
  }
  delete[] record;
  //  for(int i=0;i<splitNum;i++){
  //    close(smallFds[i]);
  //  }
  //  delete[] smallFds;
  //  H5Tclose(datatype);
  H5Dclose(datasetID);
  H5Sclose(dspaceID);
  //  H5Sclose(memspace2);
  H5Fclose(fileID);
  return nullptr;
}
