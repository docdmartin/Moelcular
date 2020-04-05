#ifndef _WRITE_VTU__
#define _WRITE_VTU__

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <stdint.h>

/*
           Cell Types
    1) vertex         - stand-alone point
    2) poly vertex    - collection of stand-alone points (no link)
    3) line           - link between two vertices
    4) poly line      - piece-wise line between ordered vertices
    5) triangle       - triangle surface using three vertices
    6) triangle strip - triangles constructed between triplets of ordered vertices

       1 - 3 - 5 - ... - N
       | \ | \ | \     \ |
       0 - 4 - 6 - ... - N-1

    7) polygon - order link between vertices, closed on first-last
    8) pixel   - square block of four vertices  (origin, X, Y, XY)
    9) quad    - just like pixel, but any arrangement (origin, X, XY, Y) - counter-clockwise
    10) tetra  - four point tetrahedron - solid
    11) voxel  - square brick using eight points (origin, X, Y, XY, Z, XZ, YZ, XYZ)
    12) hexahedron - 8 point volumn like voxel, (origin, X, XY, Y, Z, XZ, XYZ, YZ)
    13) wedge  - six point prysm (lower triangle then upper triangle)
    14) pyramid - five point solid, base first and ordered like quad


    Default options:
      mCells.empty()         -> 1) vertex, default connections are created locally
      mCells[0].size() == 1  -> 1) vertex, dummy connections provided
      mCells[0].size() == 2  -> 3) line
      mCells[0].size() == 3  -> 5) triangle
      mCells[0].size() == 4  -> 10) tetra
      mCells[0].size() == 5  -> 14) pyramid
      mCells[0].size() == 6  -> 13) wedge
      mCells[0].size() == 7  -> error
      mCells[0].size() == 8  -> 12) hexahedron, note default order is different than voxel
*/

using namespace std;

class WriteVTU{
public:
  enum VTUDataType{
    CELL,
    POINT
  };

  WriteVTU(vector< vector<double> >& points, vector< vector<int> >& connections)
  {
    mPoints = points;
    mCells  = connections;

    mNumCells  = static_cast<int>( mCells.size() );
    mNumPoints = static_cast<int>( mPoints.size() );

cout << "Number of points: " << mNumPoints << endl;
cout << "Number of connections: " << mNumCells << endl;

    if( mNumCells == 0 ) {
      mNumCells   = mNumPoints;
      mCellType = 1;
      return;
    }

    mVertPerCell = static_cast<int>( mCells[0].size() );
cout << "Vertices per connection: " << mVertPerCell << endl;


    switch( mVertPerCell ) {
      case 0:
        mVertPerCell = 1;
      case 1:
        mNumCells = mNumPoints;
        mCellType = 1;
        return;
      case 2:
        mCellType = 3;
        return;
      case 3:
        mCellType = 5;
        return;
      case 4:
        mCellType = 10;
        return;
      case 5:
        mCellType = 14;
        return;
      case 6:
        mCellType = 13;
        return;
      case 8:
        mCellType = 12;
        return;
      default:
        mCellType = 0;
        return;
    }

  }
  WriteVTU(vector< vector<double> >& points)
  {
    mPoints = points;
    mCells.clear();

    mCellType    = 1;
    mNumPoints   = static_cast<int>( mPoints.size() );
    mNumCells    = mNumPoints;
    mVertPerCell = 1;
  }
  ~WriteVTU(){}

  void AddAttributeInteger( VTUDataType dt, string attr_name, vector<int>& data )
  {
    if( dt == CELL ){
      if( data.size() != mCells.size() ) {
        cout << "Cell data provided was wrong size" << endl;
        return;
      }

      mCellMetricInt.insert( pair<string, vector<int>>(attr_name, data) );
    }
    else if( dt == POINT ){
      if( data.size() != mPoints.size() ) {
        cout << "Point data provided was wrong size" << endl;
        return;
      }

      mPointMetricInt.insert( pair<string, vector<int>>(attr_name, data) );
    }
    else {
      cout << "Data type provide was not recognized" << endl;
    }
  }
  void AddAttributeDouble ( VTUDataType dt, string attr_name, vector<double>& data ){
    if( dt == CELL ){
      if( data.size() != mCells.size() ) {
        cout << "Cell data provided was wrong size" << endl;
        return;
      }

      mCellMetricDouble.insert( pair<string, vector<double>>(attr_name, data) );
    }
    else if( dt == POINT ){
      if( data.size() != mPoints.size() ) {
        cout << "Point data provided was wrong size" << endl;
        return;
      }

      mPointMetricDouble.insert( pair<string, vector<double>>(attr_name, data) );
    }
    else {
      cout << "Data type provide was not recognized" << endl;
    }
  }

  void AddVectorField( string attr_name, vector<vector<double>>& data ) {
    if( data.size() != mPoints.size() ) {
      cout << "Vector field must be same size as point array" << endl;
      return;
    }
    mVectorField.insert( pair<string, vector<vector<double>>>( attr_name, data ) );
  }




  void WriteVTUFile( string full_filename ){
    if( mCellType == 0 ) {
      cout << "Invalid cell type" << endl;
      return;
    }

    fstream vtk_file;
    vtk_file = fstream(full_filename, ios::out | ios::binary);

    mOffset = 0;
    writeHeader( vtk_file );

    writePoints(  vtk_file );
    writeCells(   vtk_file );
    writeMetrics( vtk_file );

    writeFooter( vtk_file );
  }

private:

  vector< vector<double> > mPoints;
  vector< vector<int   > > mCells;

  uint8_t mCellType;
  int     mVertPerCell;

  map<string, vector<double>> mPointMetricDouble;
  map<string, vector<int   >> mPointMetricInt;

  map<string, vector<double>> mCellMetricDouble;
  map<string, vector<int   >> mCellMetricInt;

  map<string, vector<vector<double>>> mVectorField;

  int mOffset;
  int mNumPoints;
  int mNumCells;


  void writeHeader( fstream& vtk_file ) {

    string header_xml_version = "<?xml version=\"1.0\"?>\n";
    string header_vtk_type    = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    string start_unstr_grid   = "<UnstructuredGrid>\n";

    string define_pieces      = "<Piece NumberOfPoints=\"" + to_string(mNumPoints) + "\" NumberOfCells=\"" + to_string(mNumCells) + "\">\n";
    vtk_file.write(header_xml_version.c_str(), header_xml_version.size());
    vtk_file.write(header_vtk_type   .c_str(), header_vtk_type   .size());
    vtk_file.write(start_unstr_grid  .c_str(), start_unstr_grid  .size());
    vtk_file.write(define_pieces     .c_str(), define_pieces     .size());

    string start_points       = "<Points>\n";
    string define_points      = "<DataArray Name=\"Points\" NumberOfComponents=\"3\" type=\"Float64\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
    mOffset += 8 + mNumPoints*3*8;
    string end_points         = "</Points>\n";
    vtk_file.write(start_points      .c_str(), start_points      .size());
    vtk_file.write(define_points     .c_str(), define_points     .size());
    vtk_file.write(end_points        .c_str(), end_points        .size());

    // Required
    string start_cells        = "<Cells>\n";
    string define_connect     = "<DataArray Name=\"connectivity\" NumberOfComponents=\"1\" type=\"Int32\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
    mOffset += 8 + mNumCells*mVertPerCell*4;
    string define_offset      = "<DataArray Name=\"offsets\" NumberOfComponents=\"1\" type=\"UInt32\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
    mOffset += 8 + mNumCells*4;
    string define_num_comp    = "<DataArray Name=\"types\" NumberOfComponents=\"1\" type=\"UInt8\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
    mOffset += 8 + mNumCells;
    string end_cells          = "</Cells>\n";
    vtk_file.write(start_cells       .c_str(), start_cells       .size());
    vtk_file.write(define_connect    .c_str(), define_connect    .size());
    vtk_file.write(define_offset     .c_str(), define_offset     .size());
    vtk_file.write(define_num_comp   .c_str(), define_num_comp   .size());
    vtk_file.write(end_cells         .c_str(), end_cells         .size());


    // Assign metric to points
    map<string, vector<int   >>::iterator vit_int    = mPointMetricInt   .begin();
    map<string, vector<double>>::iterator vit_double = mPointMetricDouble.begin();
    if( !mPointMetricInt.empty() || !mPointMetricDouble.empty() ) {
      string start_block = "<PointData scalars=\"vertex_metric\">\n";
      vtk_file.write(start_block.c_str(), start_block.size());

      for(; vit_int != mPointMetricInt.end(); ++vit_int){
        string block_str = "<DataArray Name=\"" + vit_int->first + "\" NumberOfComponents=\"1\" type=\"Int32\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
        mOffset += 8 + mNumPoints*4;
        vtk_file.write(block_str.c_str(), block_str.size());
      }

      for(; vit_double != mPointMetricDouble.end(); ++vit_double){
        string block_str          = "<DataArray Name=\"" + vit_double->first + "\" NumberOfComponents=\"1\" type=\"Float64\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
        mOffset += 8 + mNumPoints*8;
        vtk_file.write(block_str.c_str(), block_str.size());
      }

      string end_block = "</PointData>\n";
      vtk_file.write(end_block.c_str(), end_block.size());
    }

    if( !mVectorField.empty() ) {
      string start_block = "<PointData vectors=\"vector_field\">\n";
      vtk_file.write(start_block.c_str(), start_block.size());

      map<string, vector<vector<double>>>::iterator vf_it = mVectorField.begin();
      for(; vf_it != mVectorField.end(); ++vf_it){
        string block_str = "<DataArray Name=\"" + vf_it->first + "\" NumberOfComponents=\"3\" type=\"Float64\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
        mOffset += 8 + mNumPoints*3*8;
        vtk_file.write(block_str.c_str(), block_str.size());
      }

      string end_block = "</PointData>\n";
      vtk_file.write(end_block.c_str(), end_block.size());
    }


    // Assign metric to elements
    map<string, vector<int   >>::iterator eit_int    = mCellMetricInt   .begin();
    map<string, vector<double>>::iterator eit_double = mCellMetricDouble.begin();
    if( !mCellMetricInt.empty() || !mCellMetricDouble.empty() ) {
      string start_block = "<CellData scalars=\"element_metric\">\n";
      vtk_file.write(start_block.c_str(), start_block.size());

      for(; eit_int != mCellMetricInt.end(); ++eit_int){
        string block_str = "<DataArray Name=\"" + eit_int->first + "\" NumberOfComponents=\"1\" type=\"Int32\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
        mOffset += 8 + mNumCells*4;
        vtk_file.write(block_str.c_str(), block_str.size());
      }

      for(; eit_double != mCellMetricDouble.end(); ++eit_double){
        string block_str          = "<DataArray Name=\"" + eit_double->first + "\" NumberOfComponents=\"1\" type=\"Float64\" format=\"appended\" offset=\"" + to_string(mOffset) + "\"/>\n";
        mOffset += 8 + mNumCells*8;
        vtk_file.write(block_str.c_str(), block_str.size());
      }

      string end_block = "</CellData>\n";
      vtk_file.write(end_block.c_str(), end_block.size());
    }


    // Required
    string end_piece          = "</Piece>\n";
    string end_unstr_grid     = "</UnstructuredGrid>\n";
    string encoding_str       = "<AppendedData encoding=\"raw\">\n_";
    vtk_file.write(end_piece         .c_str(), end_piece         .size());
    vtk_file.write(end_unstr_grid    .c_str(), end_unstr_grid    .size());
    vtk_file.write(encoding_str      .c_str(), encoding_str      .size());
  }

  void writeFooter( fstream& vtk_file ){
    string end_of_file = "\n</AppendedData>\n</VTKFile>\n";
    vtk_file.write(end_of_file.c_str(), end_of_file.size());
    vtk_file.close();
  }



  void writePoints( fstream& vtk_file ){
    int64_t char_per_block = static_cast<int64_t>(mNumPoints)*(3*8);
    vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
    for( int cnt = 0; cnt < mNumPoints; ++cnt ) {
      double x = mPoints[cnt][0];
      double y = mPoints[cnt][1];
      double z = mPoints[cnt][2];

      vtk_file.write(reinterpret_cast<char*>(&x), sizeof(double));
      vtk_file.write(reinterpret_cast<char*>(&y), sizeof(double));
      vtk_file.write(reinterpret_cast<char*>(&z), sizeof(double));
    }
  }

  void writeCells( fstream& vtk_file ){

    // define an element as its three vertices
    int64_t char_per_block = static_cast<int64_t>(mNumCells)*(mVertPerCell*4);
    vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
    for(int cnt = 0; cnt < mNumCells; ++cnt) {
      if( mVertPerCell == 1 ) {
        vtk_file.write(reinterpret_cast<char*>(&cnt), sizeof(int));
        continue;
      }
      for(int v_cnt = 0; v_cnt < mVertPerCell; ++v_cnt) {
        int i = mCells[cnt][v_cnt];
        vtk_file.write(reinterpret_cast<char*>(&i), sizeof(int));
      }
    }

    char_per_block    = static_cast<int64_t>(mNumCells*4);
    vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
    unsigned int c = 0;
    for(int cnt = 0; cnt < mNumCells; ++cnt){
      vtk_file.write(reinterpret_cast<char*>(&c), sizeof(int));
      c += mNumCells;
    }

    // Indicate element type
    char_per_block    = static_cast<int64_t>(mNumCells);
    vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
    for(int cnt = 0; cnt < mNumCells; ++cnt){
      vtk_file.write(reinterpret_cast<char*>(&mCellType), sizeof(char));
    }
  }

  void writeMetrics( fstream& vtk_file ){

    int64_t char_per_block;
    map<string, vector<int   >>::iterator vit_int;
    map<string, vector<double>>::iterator vit_double;
    map<string, vector<int   >>::iterator eit_int;
    map<string, vector<double>>::iterator eit_double;

    vit_int = mPointMetricInt.begin();
    for(; vit_int != mPointMetricInt.end(); ++vit_int){
      char_per_block = static_cast<int64_t>(mNumPoints*4);
      vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
      for( int d : vit_int->second ){
        vtk_file.write(reinterpret_cast<char*>(&d), sizeof(int));
      }
    }

    vit_double = mPointMetricDouble.begin();
    for(; vit_double != mPointMetricDouble.end(); ++vit_double){
      char_per_block = static_cast<int64_t>(mNumPoints*8);
      vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
      for( double d : vit_double->second ){
        vtk_file.write(reinterpret_cast<char*>(&d), sizeof(double));
      }
    }



    map<string, vector<vector<double>>>::iterator vf_it = mVectorField.begin();
    for(; vf_it != mVectorField.end(); ++vf_it) {
      char_per_block = static_cast<int64_t>(mNumPoints*8*3);
      vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
      for(vector<double>& vf : vf_it->second) {
        double x = vf[0];
        double y = vf[1];
        double z = vf[2];
        vtk_file.write(reinterpret_cast<char*>(&x), sizeof(double));
        vtk_file.write(reinterpret_cast<char*>(&y), sizeof(double));
        vtk_file.write(reinterpret_cast<char*>(&z), sizeof(double));
      }
    }




    eit_int = mCellMetricInt.begin();
    for(; eit_int != mCellMetricInt.end(); ++eit_int){
      char_per_block = static_cast<int64_t>(mNumCells*4);
      vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
      for( int d : eit_int->second ){
        vtk_file.write(reinterpret_cast<char*>(&d), sizeof(int));
      }
    }

    eit_double = mCellMetricDouble.begin();
    for(; eit_double != mCellMetricDouble.end(); ++eit_double){
      char_per_block = static_cast<int64_t>(mNumCells*8);
      vtk_file.write(reinterpret_cast<char*>(&char_per_block), sizeof(double));
      for( double d : eit_double->second ){
        vtk_file.write(reinterpret_cast<char*>(&d), sizeof(double));
      }
    }
  }

};

#endif
