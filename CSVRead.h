#ifndef ____CSVREAD_
#define ____CSVREAD_

#include<fstream>
#include<string>
#include<vector>
#include<map>

using namespace std;

class CSVRead {
public:
	CSVRead();
	~CSVRead();

	bool           OpenCSVFile(string);
	bool           ReadCSVHeader(vector<string>);
	vector<string> ReadCSVRecord();

private:
	vector<string> getCsvFields(string);

	ifstream       mCSVFile;
	map<int, int>  mColumnIndexMap;
};

#endif