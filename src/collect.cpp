#include<iostream>
#include<stdio.h>
#include<string>
#include<string.h>
#include<unordered_map>
#include<set>
#include<algorithm>
#include<vector>
#include<stack>
#include<cmath>
#include<sstream>
#include<deque>
using namespace std;
//************** debug information ***************
#include<fstream>
#include<dirent.h>
ofstream ofile;
time_t ntm;
time_t prtm;
//************************************************

const bool useAll = 1;
const bool notUseAll = 0;

struct tm* current_time;

typedef struct mapinfo{
	int chrm;
	int pos;
	int mapq;
	bool brevs;
	mapinfo(){pos = -1;};
	mapinfo(int c, int p, bool b, int m):chrm(c), pos(p), brevs(b), mapq(m){};
}mapinfo;

vector<vector<mapinfo> > vvmapinfo;
vector<vector<mapinfo> > vvmapinfo2;


struct Lib
{
	char FilePe[18];
	char FileSe[18];
	int MinSpan;//min_insertsize
	int MaxSpan;//max_insertsize
	int aveSpan;//ave_insertsize
	int MinOverLap;//K值
	int ReadLength;//读数长度
	int errolp;
	int std;//standard deviation
};


string Reverse_Compliment(const string& ori)
{
	string ans=ori;
	int len=ori.size();
	for(int i=len-1;i>=0;i--)
	{
		if(ori[i]=='A') ans[len-i-1]='T';
		else if(ori[i]=='C') ans[len-i-1]='G';
		else if(ori[i]=='G') ans[len-i-1]='C';
		else if(ori[i]=='T') ans[len-i-1]='A';
	}
	return ans;
}


int main(){
	//*************** code add by lzx *****************
	DIR *dir=NULL;
	dir = opendir("log");
	if(!dir){
		system("mkdir log");
	}
	
	ofile.open("log/collect.log");
	
	prtm = time(NULL); 
	ntm = prtm;
	struct tm* current_time = localtime(&ntm);
	ofile<<"program start...\t"<< ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	//************************************************
	
	int readsnum;
	ifstream ifnum("mapping/reads_num.info");
	ifnum>>readsnum;
	ifnum.close();
	vvmapinfo.resize(readsnum);
	vvmapinfo2.resize(readsnum);
	
	int readsmapped1(0), totalreads1(0);
	ifstream ifsam("mapping/single.2.1.1.sam");
	
	int chrm, seqno, flag, pos, mapq;
	int mul1(0);
	bool brevs;
	const int szbuff = 1024;
	char buff[szbuff];
	
	while(ifsam.getline(buff,szbuff)){
		stringstream ss(buff);
		ss>>seqno>>flag>>chrm>>pos>>mapq;
		//cout<<seqno<<'\t'<<flag<<'\t'<<chrm<<'\t'<<pos<<endl;
		totalreads1++;
		if(! (flag & 0x04) ){
			brevs = flag & 0x10;
			readsmapped1++;
			vvmapinfo[seqno].push_back(mapinfo(chrm, pos, brevs, mapq));
		}
	}
	ifsam.close();
	cout<<"#1 totalreads1 = "<<totalreads1<<endl;
	cout<<"#1 readsmapped1 = "<<readsmapped1<<endl;
	
	int readsmapped2(0), totalreads2(0);
	int mul2(0);
	ifsam.open("mapping/single.2.1.2.sam");
	
	while(ifsam.getline(buff,szbuff)){
		stringstream ss(buff);
		ss>>seqno>>flag>>chrm>>pos>>mapq;
		totalreads2++;
		if(! (flag & 0x04) ){
			brevs = flag & 0x10;
			readsmapped2++;
			vvmapinfo2[seqno].push_back(mapinfo(chrm, pos, brevs, mapq));
		}
	}
	ifsam.close();
	cout<<"#2 totalreads2 = "<<totalreads2<<endl;
	cout<<"#2 readsmapped2 = "<<readsmapped2<<endl;
	
	
	ofstream ofclct("mapping/collect2");
	for(int i = 0; i < vvmapinfo.size(); i++){
		if(vvmapinfo[i].size()>1){
			//cout<<vvmapinfo[i].size()<<'\t';
			// cout<<"---------"<<endl;
			// for(int p = 0; p < vvmapinfo[i].size(); p++){
				// cout<<vvmapinfo[i][p].chrm<<'\t'<<"2\t"<<vvmapinfo[i][p].pos<<'\t'<<vvmapinfo[i][p].brevs<<'\t'<<vvmapinfo[i][p].mapq<<endl;
			// }
			// cout<<"---------"<<endl;
			mul1++;
		}
		if(vvmapinfo2[i].size()>1){
			//cout<<vvmapinfo2[i].size()<<'\t';
			// cout<<"---------"<<endl;
			// for(int p = 0; p < vvmapinfo2[i].size(); p++){
				// cout<<vvmapinfo2[i][p].mapq<<endl;
			// }
			// cout<<"---------"<<endl;
			mul2++;
		}
		if(vvmapinfo[i].empty() || vvmapinfo2[i].empty())continue;
		for(int p = 0; p < vvmapinfo[i].size(); p++){
			for(int q = 0; q < vvmapinfo2[i].size(); q++){
				if(vvmapinfo[i][p].chrm == vvmapinfo2[i][q].chrm)continue;
				ofclct<<vvmapinfo[i][p].chrm<<'\t'<<"2\t"<<vvmapinfo[i][p].pos<<'\t'<<vvmapinfo[i][p].brevs<<endl;
				ofclct<<vvmapinfo2[i][q].chrm<<'\t'<<"2\t"<<vvmapinfo2[i][q].pos<<'\t'<<vvmapinfo2[i][q].brevs<<endl;
			}
		}
		
		// for(int p = vvmapinfo[i].size()-1; p < vvmapinfo[i].size(); p++){
			// for(int q = vvmapinfo[i].size()-1; q < vvmapinfo2[i].size(); q++){
				// if(vvmapinfo[i][p].chrm == vvmapinfo2[i][q].chrm)continue;
				// ofclct<<vvmapinfo[i][p].chrm<<'\t'<<"2\t"<<vvmapinfo[i][p].pos<<'\t'<<vvmapinfo[i][p].brevs<<endl;
				// ofclct<<vvmapinfo2[i][q].chrm<<'\t'<<"2\t"<<vvmapinfo2[i][q].pos<<'\t'<<vvmapinfo2[i][q].brevs<<endl;
			// }
		// }
	}
	
	cout<<endl;
	cout<<"#2 mul1 = "<<mul1<<endl;
	cout<<"#2 mul2 = "<<mul2<<endl;
	return 0;
}