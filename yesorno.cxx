#include <iostream>
#include <string>

using namespace std;

bool yesorno(const string& question){
	char answer = 'w';
	cout << question << " (y/n): ";
	cin >> answer;
	cout<< endl;
	while(answer!='y' && answer!='n' && answer!='Y' && answer!='N'){
		cout << "Answer 'y' or 'n': ";
		cin >> answer;
		cout<< endl;
	}
	if(answer=='y' || answer=='Y'){
		return true;
	}else{
		return false;
	}
};
