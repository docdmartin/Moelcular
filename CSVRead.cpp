#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

#include "CSVRead.h"

using namespace std;

CSVRead::CSVRead() {
}

CSVRead::~CSVRead() {
}

bool CSVRead::OpenCSVFile(string filename){

	mCSVFile.open(filename.c_str());
	return mCSVFile.is_open();
}

bool CSVRead::ReadCSVHeader(vector<string> header_str_vec){

	mColumnIndexMap.clear();
    
    string header_str;
    
    getline( mCSVFile, header_str );

	stringstream ss(header_str);
	vector<string> header;

	string temp_field;
	while( getline( ss, temp_field, ',' ) ){

		while(temp_field.size() > 0 && temp_field[temp_field.size()-1] < 33){
			temp_field = temp_field.substr(0, temp_field.size()-1);
		}

		if(temp_field[0] == '"' && temp_field[temp_field.size()-1] == '"'){
			temp_field = temp_field.substr(1, temp_field.size()-2);
		}
		header.push_back(temp_field);
	}

	for(int cnt = 0; cnt < static_cast<int>(header.size()); cnt++){

		for(int h = 0; h < static_cast<int>(header_str_vec.size()); h++){
			if(header[cnt].compare(header_str_vec[h]) == 0){
				mColumnIndexMap.insert(pair<int,int>(cnt,h));

				if(mColumnIndexMap.size() == header_str_vec.size()){
					return true;
				}

				break;
			}
		}
	}

    mCSVFile.close();

	return false;
}

vector<string> CSVRead::ReadCSVRecord(){

	vector<string> return_str;
    
    string record_str;
    if(getline(mCSVFile, record_str)){
        
		if(count(record_str.begin(), record_str.end(), '"')%2 == 1){
            string tmp_str;
			while(getline( mCSVFile, tmp_str )){
				record_str = record_str + tmp_str;
				if(count(record_str.begin(), record_str.end(), '"')%2 == 0)
					break;
			}
		}
        
		return_str = getCsvFields(record_str);
        
		if(return_str.size() == 0){
			mCSVFile.close();
			return return_str;
		}
        
		return return_str;
    }
    
	mCSVFile.close();
	return_str.clear();
	return return_str;
}


vector<string> CSVRead::getCsvFields(string text_string){
	stringstream ss(text_string);

	map<int, int>::iterator it;
	int curr_column = -1;

	vector<string> fields;
	fields.resize(mColumnIndexMap.size());

	int tgt_col;
	int tgt_vec;
	string temp_field;
	for(it = mColumnIndexMap.begin(); it != mColumnIndexMap.end(); ++it){
		tgt_col = it->first;
		tgt_vec = it->second;

		while(curr_column < tgt_col){
			++curr_column;
			if(!getline( ss, temp_field, ',' )){
				// Error record does not have required fields
				fields.clear();
				return fields;
			}
		}

		while(temp_field.size() > 0 && temp_field[temp_field.size()-1] < 33){
			temp_field = temp_field.substr(0, temp_field.size()-1);
		}

		if(temp_field[0] == '"' && temp_field[temp_field.size()-1] == '"'){
			temp_field = temp_field.substr(1, temp_field.size()-2);
		}
		fields[tgt_vec] = temp_field;
	}

	return fields;
}

