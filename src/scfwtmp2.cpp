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
struct spt{
	int posSelf;
	int posOther;
	int partdif;
	bool brvs;
	spt(){};
	spt(int ps,int po, int pd, int br):posSelf(ps), posOther(po), partdif(pd), brvs(br){}
};

struct target{
	int ctgno;
	int partdif;
	int gapdis;
	bool brevs;
	target(){};
	target(int c, int p, bool b):ctgno(c), partdif(p), brevs(b){};
};

typedef struct{
	vector<spt> vspt0;
	vector<spt> vspt0_1;
	vector<spt> vspt1;
	vector<spt> vspt1_1;
}dbvspt;

vector<bool> bused;
vector<int> finalctg;
vector<unordered_map<int, dbvspt> > vmpsptl;
vector<unordered_map<int, dbvspt> > vmpsptr;

vector<deque<int> > scfs;
vector<deque<int> > scfdirs;

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

int FileNum;
const int MaxFileNum=10;
Lib lib[MaxFileNum];//存储文库信息
struct Pair_Str 
{
	string oristr;
	string revstr;
};
vector<Pair_Str> OriVec;
vector<int> vlen[2];

bool sptCmp(spt e1,spt e2)
{
	return e1.posSelf > e2.posSelf;
}

bool sptCmp2(spt e1,spt e2)
{
	return e1.posOther > e2.posOther;
}

bool targetCmp(target e1,target e2)
{
	return e1.partdif > e2.partdif;
}

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

bool verifyPartner(vector<spt> &vspt, target &t, int workctg){
	int &avg = t.partdif;
	int &gapdis = t.gapdis;
	if(vspt.size()<100)return false;
	
	if(workctg == 32 && t.ctgno == 36 && vspt.size() == 2325){
		cout<<vspt.size()<<endl;
	}
	
	int sum = 0;
	//const int decrement = 200;
	const int decrement = 400;
	int stage = (vspt[0].posSelf/100)*100;
	//if(vspt[0].posSelf/100%2 == 1)stage -= 100;
	
	bool bres = true;
	for(int i = 0;i < vspt.size(); i++){
		sum += vspt[i].partdif;
		if(stage == 0){
			bres = true;
			continue;
		}
		int dif = vspt[i].posSelf - stage;
		if(dif >=0 && dif < decrement){
			stage -= decrement;
		};
		if(dif < 0){
			//cout<<vspt[i].posSelf<<'\t'<<vspt[i].partdif<<endl;
			return false;
		}
	}
	avg = sum/vspt.size();
	int qrtsum = 0;
	for(int i = 0;i < vspt.size(); i++){
		qrtsum += (avg - vspt[i].partdif)*(avg - vspt[i].partdif);
	}
	int std = sqrt(qrtsum/vspt.size());
	sum = 0;
	int numcnt = 0;
	for(int i = 0;i < vspt.size(); i++){
		if(vspt[i].partdif >= avg - 1.5*std&&vspt[i].partdif <= avg + 1.5*std)
		{
			sum += vspt[i].partdif;
			numcnt++;
		}
	}
	if(numcnt){
		avg = sum/numcnt;
		gapdis = lib[FileNum].aveSpan - avg;
	}
	int zonemin, zonemax;
	bool binzone(false);
	zonemin = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std - OriVec[t.ctgno].oristr.size();
	if(zonemin < 0)zonemin = 0;
	zonemax = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(zonemax > OriVec[workctg].oristr.size())zonemax = OriVec[workctg].oristr.size();
	zonemin = zonemin + (zonemax - zonemin) * 0.25;
	zonemax = zonemax - (zonemax - zonemin) * 0.25;
	
	cout<<"pair : "<<workctg<<'\t'<<t.ctgno<<'\t'<<gapdis<<endl;
	cout<<OriVec[workctg].oristr.size()<<'\t'<<OriVec[t.ctgno].oristr.size()<<endl;
	cout<<zonemin<<'\t'<<zonemax<<endl;
	
	if(OriVec[workctg].oristr.size()<500)binzone = true;
	if(vspt[vspt.size()-1].posSelf <= zonemin && vspt[0].posSelf >= zonemax)binzone = true;
	
	if(!binzone)return false;
	
	// int zonemin2, zonemax2;
	// bool binzone2(false);
	// zonemin2 = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std - OriVec[workctg].oristr.size();
	// if(zonemin2 < 0)zonemin2 = 0;
	// zonemax2 = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	// if(zonemax2 > OriVec[t.ctgno].oristr.size())zonemax2 = OriVec[t.ctgno].oristr.size() ;
	// zonemin2 = zonemin2 + (zonemax2 - zonemin2) * 0.25;
	// zonemax2 = zonemax2 - (zonemax2 - zonemin2) * 0.25;
	
	// cout<<zonemin2<<'\t'<<zonemax2<<endl;
	
	// if(OriVec[t.ctgno].oristr.size()<500)binzone2 = true;
	// if(vspt[vspt.size()-1].posSelf <= zonemin2 && vspt[0].posSelf >= zonemax2)binzone2 = true;
	
	//if(!binzone2)return false;
	
	bool bleft(false), bright(false);
	int minpos1, maxpos1;
	minpos1 = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(minpos1 > OriVec[workctg].oristr.size()){
		return false;
	};
	maxpos1 = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std;
	if(maxpos1 > OriVec[workctg].oristr.size()){
		maxpos1 = OriVec[workctg].oristr.size();
	};
	for(int i = 0;i < vspt.size(); i++){
		if(!bleft){
			if(vspt[i].posSelf >= minpos1 && vspt[i].posSelf <= maxpos1){
				bleft = true;
				break;
			}
		}
	}
	
	int minpos2, maxpos2;
	minpos2 = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(minpos2 > OriVec[t.ctgno].oristr.size()){
		bright = 1;
	}else{
		maxpos2 = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std;
		if(maxpos2 > OriVec[t.ctgno].oristr.size()){
			maxpos2 = OriVec[t.ctgno].oristr.size();
		};
		for(int i = vspt.size() - 1;i >= 0; i--){
			if(!bright){
				if(vspt[i].posOther >= minpos2 && vspt[i].posOther <= maxpos2){
					bright = true;
					break;
				}
			}
		}
	};
	
	if(!bleft || !bright)return false;
	
	return bres;
}

bool verify(vector<spt> &vspt, vector<spt> &vspt2, target &t, int workctg, bool bvr){
	int &avg = t.partdif;
	int &gapdis = t.gapdis;
	if(vspt.size()<50)return false;
	
	ofile<<"pair : "<<workctg<<'\t'<<t.ctgno<<'\t'<<vspt.size()<<endl;
	if(workctg == 15 && t.ctgno == 144){
		ofile<<vspt.size()<<endl;
	}
	
	int sum = 0;
	//const int decrement = 200;
	const int decrement = 400;
	int stage = (vspt[0].posSelf/100)*100;
	if(vspt[0].posSelf/100%2 == 1)stage -= 100;
	
	bool bres = true;
	for(int i = 0;i < vspt.size(); i++){
		sum += vspt[i].partdif;
		if(stage <= 0){
			bres = true;
			continue;
		}
		int dif = vspt[i].posSelf - stage;
		if(dif >=0 && dif < decrement){
			stage -= decrement;
		};
		if(dif < 0){
			ofile<<"dif < 0"<<endl;
			ofile<<vspt[i].posSelf<<'\t'<<vspt[i].partdif<<endl;
			return false;
		}
	}
	avg = sum/vspt.size();
	int qrtsum = 0;
	for(int i = 0;i < vspt.size(); i++){
		qrtsum += (avg - vspt[i].partdif)*(avg - vspt[i].partdif);
	}
	int std = sqrt(qrtsum/vspt.size());
	sum = 0;
	int numcnt = 0;
	for(int i = 0;i < vspt.size(); i++){
		if(vspt[i].partdif >= avg - 1.5*std&&vspt[i].partdif <= avg + 1.5*std)
		{
			sum += vspt[i].partdif;
			numcnt++;
		}
	}
	if(numcnt){
		avg = sum/numcnt;
		gapdis = lib[FileNum].aveSpan - avg;
	}
	int zonemin, zonemax;
	bool binzone(false);
	zonemin = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std - OriVec[t.ctgno].oristr.size();
	if(zonemin < 0)zonemin = 0;
	zonemax = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(zonemax > OriVec[workctg].oristr.size())zonemax = OriVec[workctg].oristr.size();
	zonemin = zonemin + (zonemax - zonemin) * 0.25;
	zonemax = zonemax - (zonemax - zonemin) * 0.25;
	
	ofile<<t.brevs<<'\t'<<bvr<<'\t'<<gapdis<<endl;
	ofile<<OriVec[workctg].oristr.size()<<'\t'<<OriVec[t.ctgno].oristr.size()<<endl;
	ofile<<zonemin<<'\t'<<zonemax<<endl;
	ofile<<vspt[vspt.size()-1].posSelf<<'\t'<<vspt[0].posSelf<<endl;
	
	
	if(OriVec[workctg].oristr.size()<500)binzone = true;
	if(vspt[vspt.size()-1].posSelf <= zonemin && vspt[0].posSelf >= zonemax)binzone = true;
	
	if(!binzone)return false;
	
	target tt;
	tt.ctgno = workctg;
	tt.brevs = false;
	if(t.brevs)tt.brevs = true;
	
	bool bverified;
	ofile<<"verifyPartner"<<endl;
	
	int zonemin2, zonemax2;
	bool binzone2(false);
	zonemin2 = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std - OriVec[workctg].oristr.size();
	if(zonemin2 < 0)zonemin2 = 0;
	zonemax2 = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(zonemax2 > OriVec[t.ctgno].oristr.size())zonemax2 = OriVec[t.ctgno].oristr.size() ;
	zonemin2 = zonemin2 + (zonemax2 - zonemin2) * 0.25;
	zonemax2 = zonemax2 - (zonemax2 - zonemin2) * 0.25;
	
	ofile<<zonemin2<<'\t'<<zonemax2<<endl;
	
	if(OriVec[t.ctgno].oristr.size()<500)binzone2 = true;
	ofile<<vspt2[vspt.size()-1].posOther<<'\t'<<vspt2[0].posOther<<endl;
	if(vspt2[vspt.size()-1].posOther <= zonemin2 && vspt2[0].posOther >= zonemax2)binzone2 = true;
	
	if(!binzone2){
		ofile<<"verifyPartner failed"<<endl;
		return false;
	}
	
	bool bleft(false), bright(false);
	int minpos1, maxpos1;
	minpos1 = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(minpos1 > OriVec[workctg].oristr.size()){
		return false;
	};
	maxpos1 = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std;
	if(maxpos1 > OriVec[workctg].oristr.size()){
		maxpos1 = OriVec[workctg].oristr.size();
	};
	for(int i = 0;i < vspt.size(); i++){
		if(!bleft){
			if(vspt[i].posSelf >= minpos1 && vspt[i].posSelf <= maxpos1){
				bleft = true;
				break;
			}
		}
	}
	
	int minpos2, maxpos2;
	minpos2 = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(minpos2 > OriVec[t.ctgno].oristr.size()){
		bright = 1;
	}else{
		maxpos2 = lib[FileNum].aveSpan - gapdis + 3*lib[FileNum].std;
		if(maxpos2 > OriVec[t.ctgno].oristr.size()){
			maxpos2 = OriVec[t.ctgno].oristr.size();
		};
		for(int i = vspt.size() - 1;i >= 0; i--){
			if(!bright){
				if(vspt[i].posOther >= minpos2 && vspt[i].posOther <= maxpos2){
					bright = true;
					break;
				}
			}
		}
	};
	
	if(!bleft || !bright)return false;
	
	return bres;
}

void join(int workctg, vector<target> &vtgt, int dir){
	ofile<<"joined : ";
	int prixlen = 0;
	for(int i=0; i<vtgt.size(); i++){
		ofile<<"("<<vtgt[i].ctgno<<'\t'<<vtgt[i].brevs<<")\t";;
		cout<<"("<<vtgt[i].ctgno<<'\t'<<vtgt[i].brevs<<")\t"<<endl;
		bused[vtgt[i].ctgno] = 1;
		int noOfn = (vtgt[i].gapdis - prixlen);
		prixlen += (noOfn + OriVec[vtgt[i].ctgno].oristr.size());
		if(noOfn <0 )noOfn = 0;
		if(dir == 1){//right extension
			scfs[workctg].push_back(noOfn);
			scfdirs[workctg].push_back(noOfn);
			// if(!vtgt[i].brevs){
				// for(int r = 0; r < scfs[abs(vtgt[i].ctgno)].size(); r++){
					// scfs[workctg].push_back(scfs[abs(vtgt[i].ctgno)][r]);
					// scfdirs[workctg].push_back(scfdirs[abs(vtgt[i].ctgno)][r]);
				// }
			// }else{
				// for(int r =  ((int) scfs[abs(vtgt[i].ctgno)].size()) -1; r >= 0; r--){
					// scfs[workctg].push_back(scfs[abs(vtgt[i].ctgno)][r]);
					// scfdirs[workctg].push_back((-1)scfdirs[abs(vtgt[i].ctgno)][r]);
				// }
			// }
			
			scfs[workctg].push_back(abs(vtgt[i].ctgno));
			scfdirs[workctg].push_back(vtgt[i].brevs);
		}else{
			scfs[workctg].push_front(noOfn);
			if(vtgt[i].brevs)scfs[workctg].push_front(vtgt[i].ctgno*-1);
			else scfs[workctg].push_front(vtgt[i].ctgno);
		}
	}
	ofile<<endl;
}

void collectTargets(int workctg, vector<target> &vtgt, int rightedge, bool brevs, bool buseAll){
	target t;
	int i = workctg;
	//right extension  unordered_map<int, vector<spt> >
	ofile<<"rightedge = "<<rightedge<<endl;
	
	if(!brevs){
		ofile<<"candidates  with surpports = "<<vmpsptl[rightedge].size()<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptl[rightedge].begin(); it!=vmpsptl[rightedge].end(); it++){
			if(!buseAll && bused[it->first])continue;
			ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
			t.ctgno = it->first;
			t.brevs = 0;
			if(verify(it->second.vspt0, it->second.vspt0_1, t, rightedge, false)){
				t.brevs = 0;
				vtgt.push_back(t);
			};
			
			t.brevs = 1;
			ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
			if(verify(it->second.vspt1, it->second.vspt1_1, t, rightedge, false)){
				t.brevs = 1;
				vtgt.push_back(t);
			};
			
		}
	}else{
		ofile<<"candidates  with surpports = "<<vmpsptr[rightedge].size()<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptr[rightedge].begin(); it!=vmpsptr[rightedge].end(); it++){
			if(!buseAll && bused[it->first])continue;
			ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
			t.ctgno = it->first;
			t.brevs = 1;
			if(rightedge == 36 && t.ctgno == 32){
				cout<<"here"<<endl;
			}
			if(verify(it->second.vspt0, it->second.vspt0_1, t, rightedge, true)){
				t.brevs = 1;
				vtgt.push_back(t);
			}
			ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
			t.brevs = 0;
			if(verify(it->second.vspt1, it->second.vspt1_1, t, rightedge,true )){
				t.brevs = 0;
				vtgt.push_back(t);
			}
		}
	}
	sort(vtgt.begin(), vtgt.end(), targetCmp);
}


bool scf(vector<target> &vtgt,int workctg, int dir){
	if(!vtgt.size())return false;
	int prixlen = 0;
	int rightedge;
	bool brevs;
	int validpath = 0;
	for(int i=0; i < vtgt.size()-1; i++){
		int noOfn = (vtgt[i].gapdis - prixlen);
		prixlen += (noOfn + OriVec[vtgt[i].ctgno].oristr.size());
		if(prixlen > lib[FileNum].MaxSpan){
			rightedge = vtgt[0].ctgno;
			brevs = !vtgt[0].brevs;
			vector<target> vtgt2;
			collectTargets(rightedge, vtgt2, rightedge, brevs, useAll);
			for(int i = 0; i< vtgt2.size(); i++){
				vector<target> vtgt3;
				collectTargets(vtgt2[i].ctgno, vtgt3, vtgt2[i].ctgno, !vtgt2[i].brevs, useAll);
				if(vtgt3.size() == 1 && vtgt3[0].ctgno == vtgt[0].ctgno && vtgt3[0].ctgno != workctg){
					vtgt.erase(vtgt.begin(), vtgt.begin()+1);
					return scf(vtgt, workctg, dir);
				}
			}
			return false;
		}
	}
	join(workctg, vtgt, dir);
	return true;
}

int main(){
	//*************** code add by lzx *****************
	DIR *dir=NULL;
	dir = opendir("log");
	if(!dir){
		system("mkdir log");
	}
	
	ofile.open("log/scfwtmp2.log");
	
	prtm = time(NULL); 
	ntm = prtm;
	struct tm* current_time = localtime(&ntm);
	ofile<<"program start...\t"<< ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	//************************************************
	
	if(!freopen("lib.info","r",stdin)){
		perror("lib.info");
	};
	scanf("%d",&FileNum);
	scanf("%d",&lib[0].MinOverLap);
	for(int i=1;i<=FileNum;i++)
	{
		scanf("%s%d%d%d%d%d",lib[i].FilePe,&lib[i].MinSpan,&lib[i].MaxSpan,&lib[i].std,&lib[i].ReadLength, &lib[i].errolp);
		lib[i].MinOverLap = lib[0].MinOverLap;
		lib[i].aveSpan=(lib[i].MinSpan+lib[i].MaxSpan)/2;
		lib[i].MaxSpan=lib[i].aveSpan+3*lib[i].std;
		lib[i].MinSpan=lib[i].aveSpan-3*lib[i].std;
		sprintf(lib[i].FilePe,"pe_%d.txt",i);
		sprintf(lib[i].FileSe,"se_%d.txt",i);
	}
	
	if(!freopen("combine.fasta","r",stdin)){
		perror("combine.fasta");
	};
	string str;
	Pair_Str ps;
	while(cin>>str)
	{
        cin>>str;
		ps.oristr=str;
		ps.revstr=Reverse_Compliment(ps.oristr);
		OriVec.push_back(ps);
		scfs.push_back(deque<int>());
		scfdirs.push_back(deque<int>());
		scfs[OriVec.size()-1].push_back(OriVec.size()-1);
		scfdirs[OriVec.size()-1].push_back(0);
	}
	bused.resize(OriVec.size());
	vmpsptl.resize(OriVec.size());
	vmpsptr.resize(OriVec.size());
	vlen[0].resize(OriVec.size());
	vlen[1].resize(OriVec.size());
	
	for(int i=0;i<bused.size();i++)bused[i] = false;
	
	ifstream ifbs("mapping/bases2.fa");
	while(ifbs>>str)
	{
		bool bsplit = 0;
		for(int i=0;i<str.size();i++){
			if(str[i] == '>'){
				str[i] = ' ';
			}
			if(str[i] == ':'){
				str[i] = '\t';
				bsplit = 1;
				break;
			}
		}
		stringstream ss;
		int ctgno,flk;
		if(bsplit){
			ss<<str;
			ss>>ctgno>>flk;
			ifbs>>str;
			vlen[flk][ctgno] = str.size();
		}else{
			ss<<str;
			ss>>ctgno;
			ifbs>>str;
			vlen[0][ctgno] = str.size();
			vlen[1][ctgno] = str.size();
		}
	}
	
	ifstream ifs("mapping/collect2");
	
	//seq [no] 0:left flank, 1:right flank, 2:full contig
	
	int ctgno1, flk1, pos1, dir1;
	int ctgno2, flk2, pos2, dir2;
	int flkdf1, flkdf2;
	
	int spts = 0;
	while(ifs>>ctgno1>>flk1>>pos1>>dir1){
		ifs>>ctgno2>>flk2>>pos2>>dir2;
		spts++;
		if(ctgno1 == 0 && ctgno2 == 115 || ctgno1 == 115 && ctgno2 == 0){
			// cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
			// cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
			// cout<<"-------------------------"<<endl;
			// continue;
		}
		
		if(ctgno1>OriVec.size()||ctgno2>OriVec.size())cout<<OriVec.size()<<'\t'<<ctgno1<<'\t'<<ctgno2<<endl;
		//continue;
		//lib[FileNum].MaxSpan = 100000000;
		if(dir1){
			if(dir2){ //dif strand 1 -> -2
				//if(flk1 == 0 || flk2 == 1)continue;
				flkdf1 = vlen[!!flk1][ctgno1]-pos1;
				flkdf2 = vlen[!!flk2][ctgno2]-pos2;
				if(flkdf1<0||flkdf2<0){
					cout<<"1"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno1][ctgno2].vspt1.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 1)); 
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno2][ctgno1].vspt1.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 1));
			}else{ //same strand 1 -> 2
				//if(flk1 == 0 || flk2 == 1)continue;
				flkdf1 = vlen[!!flk1][ctgno1]-pos1;
				flkdf2 = pos2+1; 
				if(flkdf1<0||flkdf2<0){
					cout<<"2"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				// if(ctgno1 == 0 && ctgno2 == 8 || ctgno1 == 8 && ctgno2 == 0){
					// cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					// cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					// cout<<flkdf1<<'\t'<<flkdf2<<endl;
					// cout<<flkdf1+flkdf2<<endl;
				// }
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno1][ctgno2].vspt0.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 0)); 
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno2][ctgno1].vspt0.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 0));
			}
		}else{
			if(dir2){ // same strand 2 -> 1
				//if(flk1 == 1 || flk2 == 0)continue;
				flkdf1 = pos1+1;
				flkdf2 = vlen[!!flk2][ctgno2]-pos2;
				if(flkdf1<0||flkdf2<0){
					cout<<"3"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno2][ctgno1].vspt0.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 0));
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno1][ctgno2].vspt0.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 0));
			}else{ //dif strand -2 -> 1
				//if(flk1 == 1 || flk2 == 0)continue;
				flkdf1 = pos1+1;
				flkdf2 = pos2+1;
				if(flkdf1<0||flkdf2<0){
					cout<<"4"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno1][ctgno2].vspt1.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 1)); 
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno2][ctgno1].vspt1.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 1));
			}
		}
	}
	//exit(0);
	cout<<"spts = "<<spts<<endl;
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"read over, "<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	prtm = ntm; 
	//************************************************
	
	int spts2 = 0;
	int gt500 = 0;
	for(int i=0;i<OriVec.size();i++){
		
		ofile<<"==================================\n==================================\nvmpsptl:"<<endl;
		ofile<<"contig "<<i<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptl[i].begin(); it!=vmpsptl[i].end(); it++){
			sort(it->second.vspt0.begin(),it->second.vspt0.end(),sptCmp);
			sort(it->second.vspt1.begin(),it->second.vspt1.end(),sptCmp);
			it->second.vspt0_1 = it->second.vspt0;
			it->second.vspt1_1 = it->second.vspt1;
			sort(it->second.vspt0_1.begin(),it->second.vspt0_1.end(),sptCmp2);
			sort(it->second.vspt1_1.begin(),it->second.vspt1_1.end(),sptCmp2);
			
			spts2 += it->second.vspt0.size();
			spts2 += it->second.vspt1.size();
			if(it->second.vspt0.size()>300)gt500++;
			if(it->second.vspt1.size()>300)gt500++;
			
			if(it->second.vspt0.size()>0){
				ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
				for(vector<spt>::iterator itv = it->second.vspt0.begin();itv != it->second.vspt0.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			
			if(it->second.vspt1.size()>0){
				ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
				for(vector<spt>::iterator itv = it->second.vspt1.begin();itv != it->second.vspt1.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			
		}
		ofile<<"\n==================================\n==================================\n"<<endl<<endl;
		
		ofile<<"==================================\n==================================\nvmpsptr:"<<endl;
		ofile<<"contig "<<i<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptr[i].begin(); it!=vmpsptr[i].end(); it++){
			sort(it->second.vspt0.begin(),it->second.vspt0.end(),sptCmp);
			sort(it->second.vspt1.begin(),it->second.vspt1.end(),sptCmp);
			it->second.vspt0_1 = it->second.vspt0;
			it->second.vspt1_1 = it->second.vspt1;
			sort(it->second.vspt0_1.begin(),it->second.vspt0_1.end(),sptCmp2);
			sort(it->second.vspt1_1.begin(),it->second.vspt1_1.end(),sptCmp2);
			
			spts2 += it->second.vspt0.size();
			spts2 += it->second.vspt1.size();
			if(it->second.vspt0.size()>300)gt500++;
			if(it->second.vspt1.size()>300)gt500++;
			
			if(it->second.vspt0.size()>0){
				ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
				for(vector<spt>::iterator itv = it->second.vspt0.begin();itv != it->second.vspt0.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			
			if(it->second.vspt1.size()>0){
				ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
				for(vector<spt>::iterator itv = it->second.vspt1.begin();itv != it->second.vspt1.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
		}
		ofile<<"==================================\n==================================\n"<<endl<<endl;
	}
	
	cout<<"spts2 = "<<spts2<<endl;
	cout<<"gt300 = "<<gt500<<endl;
	//exit(0);
	
	int nwhile = 0;
	time_t ntm2,prtm2,ntm3,prtm3;
	ntm2 = time(NULL);
	ntm3 = time(NULL);
	int lastno, no;
	while(1)
	{
		nwhile++;
		lastno = 0;
		no = 0;
		
		for(int i=0;i<OriVec.size();i++)if(!bused[i])lastno++;
		
		for(int i=0;i<OriVec.size();i++){
			if(!bused[i])
				ofile<<"(i, nContigs) = "<<"("<<i<<", "<<OriVec.size()<<")"<<endl;
			else continue;
			bused[i] = true;
			finalctg.push_back(i);
			bool bRightDone = false;
			while(1){
				vector<target> vtgt;
				
				//right extension  unordered_map<int, vector<spt> >
				int rightedge = scfs[i][scfs[i].size()-1];
				bool brevs = scfdirs[i][scfdirs[i].size()-1];
				
				collectTargets(i, vtgt, rightedge, brevs, notUseAll);
				
				ofile<<"final candidates are: \n";
				for(int j=0;j<vtgt.size();j++){
					ofile<<vtgt[j].ctgno<<'\t'<<vtgt[j].partdif<<'\t'<<vtgt[j].gapdis<<'\t'<<vtgt[j].brevs<<'\t'<<OriVec[vtgt[j].ctgno].oristr.size()<<endl;
				}
				if(!scf(vtgt, i, 1)){
					if(bRightDone){
						ofile<<"left extension over!"<<endl;
						break;
					}else{
						bRightDone = true;
						ofile<<"right extension over!"<<endl;
						reverse(scfs[i].begin(), scfs[i].end());
						reverse(scfdirs[i].begin(), scfdirs[i].end());
						for(int j = 0;j<scfdirs[i].size();j+=2){
							scfdirs[i][j] = !scfdirs[i][j];
						}
					}
				}
			}
			
			
			
			
			//************** debug information ***************
			prtm2 = ntm2; 
			ntm2 = time(NULL);
			current_time = localtime(&ntm2);
			ofile<<"(i, nContigs) = "<<"("<<i<<", "<<OriVec.size()<<"), "<<ntm2-prtm2 <<" (s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
			ofile<<"------------------------------------------"<<endl;
			//************************************************
		}
		
		for(int i=0;i<OriVec.size();i++)if(!bused[i])no++;
		//************** debug information ***************
		prtm3 = ntm3; 
		ntm3 = time(NULL);
		current_time = localtime(&ntm3);
		ofile<<"scaffold and gap round "<<nwhile<<" over,\tscaffold success "<<lastno-no<<" times, "<<ntm3-prtm3 <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
		ofile<<"------------------------------------------"<<endl;
		//************************************************
		//if(0 == lastno-no)break;
		break;
	} 
	
	ofstream ofscf("scaffold.fasta");
	ofstream oforder("order2.info");
	int cnt = 0;
	for(int q = 0;q < finalctg.size(); q++){
		int r = finalctg[q];
		ofscf<<">"<<cnt++<<endl;
		oforder<<"scafffold "<<r<<" :"<<endl;
		for(int h= 0;h < scfs[r].size();h++){
			if(h%2 == 0){ // contig
				oforder<<scfs[r][h];
				if(scfdirs[r][h]) {
					oforder<<'-';
					ofscf<<OriVec[scfs[r][h]].revstr;
				}else {
					oforder<<'+';
					ofscf<<OriVec[scfs[r][h]].oristr;
				}
				oforder<<'\t';
			}else{ // gap N
				while(scfdirs[r][h]--)ofscf<<'N';
			}
		}
		oforder<<endl;
		ofscf<<endl;
	}
	
	return 0;
}