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
	deque<spt>  *pdspt;
	target(){};
	target(int c, int p, bool b):ctgno(c), partdif(p), brevs(b){};
};

typedef struct{
	deque<spt> vspt0;
	deque<spt> vspt0_1; //和vspt0的内容一样，只是元素顺序不一样
	deque<spt> vspt1;  //目标是反向互补的
	deque<spt> vspt1_1;
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

void collectTargets(int workctg, deque<target> &vtgt, int rightedge, bool brevs, bool buseAll, bool need2verifyParterner);

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
	return e1.gapdis < e2.gapdis;
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

const int substrLen = 20;
const int limitlen = 650;
size_t isOverlap(string &str1, string &str2){
	string sub1 = str1.substr(str1.size() - substrLen);
	string sub2 = str2.substr(0, 1000);
	return sub2.find(sub1);
}

bool verifyParterner(int workctg, int rightedge, bool brevs){
	ofile<<"verifyParterner"<<endl;
	deque<target> vtgt;
	collectTargets(rightedge, vtgt, rightedge, !brevs, useAll, false);
	//if(!isPathOK(vtgt, rightedge, rightedge))return false;
	unordered_map<int, int> umpPos;
	int prixlen = 0;
	for(int i = (int) scfs[workctg].size() - 1; i >= 0; i--){
		int ctgno = scfs[workctg][i];
		umpPos[ctgno] = prixlen * (scfdirs[workctg][i]?1:-1);
		if(--i >= 0)prixlen += scfs[workctg][i];
		else break;
	}
	
	int startpos = 0;
	int maxid = 0;
	
	for(int i = 0; i < vtgt.size(); i++){
		if(vtgt[i].pdspt->size() > vtgt[maxid].pdspt->size()) maxid = i;
	}
	
	if(umpPos.count(vtgt[maxid].ctgno)){
		ofile<<"verifyParterner success"<<endl;
		return true;
	}else return false;
}

bool sizeSatisfied(int c1, int c2, int sz){ //c1 is rightedge, c2 is candidate
	//if( (OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && vspt.size()<20)return false;
	//if( !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && vspt.size()<10)return false;
	
}

int getGapdis(int start, int end, deque<spt> &vspt, target &t, int workctg, int origctg, bool bvr){
	int &avg = t.partdif;
	int gapdis;
	
	int sum = 0;
	int stage = (vspt[0].posSelf/100)*100;
	if(vspt[0].posSelf/100%2 == 1)stage -= 100;
	
	for(int i = start;i < end; i++){
		sum += vspt[i].partdif;
	}
	avg = sum/(end - start);
	int qrtsum = 0;
	for(int i = start;i < end; i++){
		qrtsum += (avg - vspt[i].partdif)*(avg - vspt[i].partdif);
	}
	int std = sqrt(qrtsum/vspt.size());
	sum = 0;
	int numcnt = 0;
	for(int i = start;i < end; i++){
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
	return gapdis;
}

bool verifyzone(int start, int end, deque<spt> &vspt, target &t, int workctg, int origctg, bool bvr, bool need2verifyParterner){
	int &avg = t.partdif;
	int &gapdis = t.gapdis;
	
	ofile<<"pair : "<<workctg<<'\t'<<t.ctgno<<'\t'<<vspt.size()<<endl;
	ofile<<"zone : "<<start<<'\t'<<end<<endl;
	
	vector<int > vgapdis;
	vector<int > vgdsSpt;
	end--;
	int i = end - 1;
	for(; i >= 0; i--){
		if(vspt[i].posSelf - vspt[i+1].posSelf > 2*lib[FileNum].std){
			//find a zone, verify this zone
			if((OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 20 || !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 10){
				vgapdis.push_back(getGapdis(i+1, end+1, vspt, t, workctg, origctg, bvr));
				vgdsSpt.push_back(end - i);
			}
			end = i;
		}
	}
	if((OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 20 || !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 10){
		vgapdis.push_back(getGapdis(i+1, end+1, vspt, t, workctg, origctg, bvr));
		vgdsSpt.push_back(end - i);
	}
	
	if(vgapdis.empty()){
		//vgapdis.push_back(getGapdis(0, vspt.size(), vspt, t, workctg, origctg, bvr));
		return false;
	}
	ofile<<"gapdis = ";
	for(int i = 0; i< vgapdis.size(); i++)ofile<<vgapdis[i]<<' ';
	
	int longestzone = 0;
	for(int p = 0; p < vgdsSpt.size(); p++){
		if(vgdsSpt[p] > vgdsSpt[longestzone])longestzone = p;
	}
	t.gapdis = vgapdis[longestzone];
	ofile<<", final gapdis = "<<vgapdis[longestzone]<<endl;
	if(need2verifyParterner && !verifyParterner(origctg, t.ctgno, t.brevs))return false;
	return true;
}

bool verifyStep1(deque<spt> &vspt, target &t, int workctg, int origctg, bool bvr, bool need2verifyParterner){
	if(workctg == 230){
		workctg = 230;
	}
	if( (OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && vspt.size()<20)return false;
	if( !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && vspt.size()<10)return false;
	ofile<<"verifyStep1 : "<<workctg<<'\t'<<bvr<<'\t'<<t.ctgno<<'\t'<<origctg<<endl;
	for(int p = 0; p < vspt.size(); p++){
		ofile<<'('<<vspt[p].posSelf<<','<<vspt[p].posOther<<','<<vspt[p].partdif<<"),";
	}
	ofile<<endl;
	if(verifyzone(0, vspt.size(), vspt, t, workctg, origctg, bvr, need2verifyParterner))return true;
	return false;
}

bool verifyzone2(int start, int end, deque<spt> &vspt, target &t, int workctg, int origctg, bool bvr, bool need2verifyParterner){
	int &avg = t.partdif;
	int &gapdis = t.gapdis;
	
	ofile<<"pair : "<<workctg<<'\t'<<t.ctgno<<'\t'<<vspt.size()<<endl;
	ofile<<"zone : "<<start<<'\t'<<end<<endl;
	
	int sum = 0;
	//const int decrement = 200;
	int stage = (vspt[0].posSelf/100)*100;
	if(vspt[0].posSelf/100%2 == 1)stage -= 100;
	
	for(int i = start;i < end; i++){
		sum += vspt[i].partdif;
	}
	avg = sum/(end - start);
	int qrtsum = 0;
	for(int i = start;i < end; i++){
		qrtsum += (avg - vspt[i].partdif)*(avg - vspt[i].partdif);
	}
	int std = sqrt(qrtsum/vspt.size());
	sum = 0;
	int numcnt = 0;
	for(int i = start;i < end; i++){
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
	int worklen = 0;
	if(workctg == 141){
		workctg = 141;
	}
	for(int i = scfs[origctg].size() - 1; i >= 0;){
		worklen += OriVec[scfs[origctg][i]].oristr.size();
		i--;
		if(worklen > lib[FileNum].aveSpan)break;
		if(i >= 0)worklen += scfs[origctg][i];
		i--;
	}
	zonemin = lib[FileNum].aveSpan - gapdis - OriVec[t.ctgno].oristr.size() + 3*lib[FileNum].std;
	if(zonemin < 0)zonemin = lib[FileNum].std;
	zonemax = lib[FileNum].aveSpan - gapdis - 3*lib[FileNum].std;
	if(zonemax > worklen)zonemax = worklen - lib[FileNum].std;
	//zonemin = zonemin + (zonemax - zonemin) * 0.25;
	//zonemax = zonemax - (zonemax - zonemin) * 0.25;
	
	ofile<<t.brevs<<'\t'<<bvr<<'\t'<<gapdis<<endl;
	ofile<<worklen<<'\t'<<OriVec[workctg].oristr.size()<<'\t'<<OriVec[t.ctgno].oristr.size()<<endl;
	ofile<<zonemin<<'\t'<<zonemax<<endl;
	ofile<<vspt[end-1].posSelf<<'\t'<<vspt[start].posSelf<<endl;
	
	
	if(worklen<500)binzone = true;
	if(vspt[end-1].posSelf <= zonemin && vspt[start].posSelf >= zonemax)binzone = true;
	
	if(!binzone){
		ofile<<"not passed"<<endl;
		return false;
	}
	ofile<<"passed"<<endl;
	if(need2verifyParterner && !verifyParterner(origctg, t.ctgno, t.brevs))return false;
	return true;
}

bool verifyStep2(deque<spt> &vspt, target &t, int workctg, int origctg, bool bvr, bool need2verifyParterner){
	if( (OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && vspt.size()<20)return false;
	if( !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && vspt.size()<10)return false;
	//divide zones
	
	int end = vspt.size() - 1;
	int i = vspt.size() - 2;
	bool binctg = true;
	int pos = scfs[origctg].size() - 1;
	int bound = OriVec[scfs[origctg][pos]].oristr.size();
	for(; i >= 0; i--){
		if(vspt[i].posSelf > bound){
			binctg = !binctg;
			if(binctg)bound += OriVec[scfs[origctg][--pos]].oristr.size();
			else bound += abs(scfs[origctg][--pos]);
		}
		if(!binctg)continue;
		if(vspt[i].posSelf - vspt[i+1].posSelf > 3*lib[FileNum].std){
			//find a zone, verify this zone
			if((OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 20 || !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 10){
				if(verifyzone2(i+1, end+1, vspt, t, workctg, origctg, bvr, need2verifyParterner))return true;
			}
			end = i;
		}
	}
	if((OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 20 || !(OriVec[t.ctgno].oristr.size() > limitlen && OriVec[workctg].oristr.size() > limitlen) && end - i >= 10){
		if(verifyzone2(i+1, end+1, vspt, t, workctg, origctg, bvr, need2verifyParterner))return true;
	}
	return false;
}

void merge(int w, bool bwr, int t, bool btr, int gapdis){
	ofile<<"merge : "<<w<<' '<<bwr<<' '<<t<<' '<<btr<<' '<<gapdis<<endl;
	if(w == 72 && t == 110){
		w = 72;
	}
	if(btr){
		if(bwr){
			for(unordered_map<int, dbvspt>::iterator it = vmpsptr[w].begin(); it != vmpsptr[w].end(); it++){
				if(it->first == t || bused[it->first])continue;
				for(int i = 0; i < it->second.vspt0.size(); i++){
					it->second.vspt0[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt0[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
				for(int i = 0; i < it->second.vspt1.size(); i++){
					it->second.vspt1[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt1[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
			}
			for(unordered_map<int, dbvspt>::iterator it = vmpsptr[w].begin(); it != vmpsptr[w].end(); it++){
				if(it->first == t || bused[it->first])continue;
				vmpsptr[t][it->first].vspt0.insert(vmpsptr[t][it->first].vspt0.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptr[t][it->first].vspt1.insert(vmpsptr[t][it->first].vspt1.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				vmpsptr[t][it->first].vspt0_1.insert(vmpsptr[t][it->first].vspt0_1.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptr[t][it->first].vspt1_1.insert(vmpsptr[t][it->first].vspt1_1.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				
				sort(vmpsptr[t][it->first].vspt0_1.begin(),vmpsptr[t][it->first].vspt0_1.end(),sptCmp2);
				sort(vmpsptr[t][it->first].vspt1_1.begin(),vmpsptr[t][it->first].vspt1_1.end(),sptCmp2);
			}
		}else{
			for(unordered_map<int, dbvspt>::iterator it = vmpsptl[w].begin(); it != vmpsptl[w].end(); it++){
				if(it->first == t)continue;
				for(int i = 0; i < it->second.vspt0.size(); i++){
					it->second.vspt0[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt0[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
				for(int i = 0; i < it->second.vspt1.size(); i++){
					it->second.vspt1[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt1[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
			}
			for(unordered_map<int, dbvspt>::iterator it = vmpsptl[w].begin(); it != vmpsptl[w].end(); it++){
				if(it->first == t)continue;
				vmpsptr[t][it->first].vspt1.insert(vmpsptr[t][it->first].vspt1.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptr[t][it->first].vspt0.insert(vmpsptr[t][it->first].vspt0.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				vmpsptr[t][it->first].vspt1_1.insert(vmpsptr[t][it->first].vspt1_1.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptr[t][it->first].vspt0_1.insert(vmpsptr[t][it->first].vspt0_1.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				
				sort(vmpsptr[t][it->first].vspt0_1.begin(),vmpsptr[t][it->first].vspt0_1.end(),sptCmp2);
				sort(vmpsptr[t][it->first].vspt1_1.begin(),vmpsptr[t][it->first].vspt1_1.end(),sptCmp2);
			}
		}
	}else{
		if(bwr){
			for(unordered_map<int, dbvspt>::iterator it = vmpsptr[w].begin(); it != vmpsptr[w].end(); it++){
				if(it->first == t)continue;
				if(it->first == 55){
					w = 110;
				}
				for(int i = 0; i < it->second.vspt0.size(); i++){
					it->second.vspt0[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt0[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
				for(int i = 0; i < it->second.vspt1.size(); i++){
					it->second.vspt1[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt1[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
			}
			for(unordered_map<int, dbvspt>::iterator it = vmpsptr[w].begin(); it != vmpsptr[w].end(); it++){
				if(it->first == t)continue;
				vmpsptl[t][it->first].vspt0.insert(vmpsptl[t][it->first].vspt0.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				vmpsptl[t][it->first].vspt1.insert(vmpsptl[t][it->first].vspt1.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptl[t][it->first].vspt0_1.insert(vmpsptl[t][it->first].vspt0_1.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				vmpsptl[t][it->first].vspt1_1.insert(vmpsptl[t][it->first].vspt1_1.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				
				sort(vmpsptl[t][it->first].vspt0_1.begin(),vmpsptl[t][it->first].vspt0_1.end(),sptCmp2);
				sort(vmpsptl[t][it->first].vspt1_1.begin(),vmpsptl[t][it->first].vspt1_1.end(),sptCmp2);
			}
		}else{
			//vmpsptl[t]
			for(unordered_map<int, dbvspt>::iterator it = vmpsptl[w].begin(); it != vmpsptl[w].end(); it++){
				if(it->first == t)continue;
				for(int i = 0; i < it->second.vspt0.size(); i++){
					it->second.vspt0[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt0[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
				for(int i = 0; i < it->second.vspt1.size(); i++){
					it->second.vspt1[i].posSelf += (gapdis + OriVec[t].oristr.size());
					it->second.vspt1[i].partdif += (gapdis + OriVec[t].oristr.size());
				}
			}
			for(unordered_map<int, dbvspt>::iterator it = vmpsptl[w].begin(); it != vmpsptl[w].end(); it++){
				if(it->first == t)continue;
				vmpsptl[t][it->first].vspt0.insert(vmpsptl[t][it->first].vspt0.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptl[t][it->first].vspt1.insert(vmpsptl[t][it->first].vspt1.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				vmpsptl[t][it->first].vspt0_1.insert(vmpsptl[t][it->first].vspt0_1.begin(), it->second.vspt0.begin(), it->second.vspt0.end());
				vmpsptl[t][it->first].vspt1_1.insert(vmpsptl[t][it->first].vspt1_1.begin(), it->second.vspt1.begin(), it->second.vspt1.end());
				
				sort(vmpsptl[t][it->first].vspt0_1.begin(),vmpsptl[t][it->first].vspt0_1.end(),sptCmp2);
				sort(vmpsptl[t][it->first].vspt1_1.begin(),vmpsptl[t][it->first].vspt1_1.end(),sptCmp2);
			}
		}
	}
}

void join(int workctg, deque<target> &vtgt, int dir){
	ofile<<"joined : ";
	int prixlen = 0;
	for(int i=0; i<vtgt.size(); i++){
		ofile<<"("<<vtgt[i].ctgno<<'\t'<<vtgt[i].brevs<<")\t";
		int noOfn = (vtgt[i].gapdis - prixlen);
		vtgt[i].gapdis = noOfn;
		prixlen += (noOfn + OriVec[vtgt[i].ctgno].oristr.size());
		if(noOfn <0 )noOfn = 0;
		if(dir == 1){//right extension
			if(noOfn < 200){
				size_t r;
				if(!scfdirs[workctg][scfdirs[workctg].size()-1]){
					if(!vtgt[i].brevs){
						r = isOverlap(OriVec[scfs[workctg][scfs[workctg].size()-1]].oristr, OriVec[abs(vtgt[i].ctgno)].oristr);
					}else{
						r = isOverlap(OriVec[scfs[workctg][scfs[workctg].size()-1]].oristr, OriVec[abs(vtgt[i].ctgno)].revstr);
					}
				}else{
					if(!vtgt[i].brevs){
						r = isOverlap(OriVec[scfs[workctg][scfs[workctg].size()-1]].revstr, OriVec[abs(vtgt[i].ctgno)].oristr);
					}else{
						r = isOverlap(OriVec[scfs[workctg][scfs[workctg].size()-1]].revstr, OriVec[abs(vtgt[i].ctgno)].revstr);
					}
				}
				if(r != string::npos){
					//cout<<scfs[workctg][scfs[workctg].size()-1]<<'\t'<<abs(vtgt[i].ctgno)<<"\tis overlapped "<<r<<endl;
				}
			}
			scfs[workctg].push_back(noOfn);
			scfdirs[workctg].push_back(noOfn);
			
			scfs[workctg].push_back(abs(vtgt[i].ctgno));
			scfdirs[workctg].push_back(vtgt[i].brevs);
		}else{
			scfs[workctg].push_front(noOfn);
			if(vtgt[i].brevs)scfs[workctg].push_front(vtgt[i].ctgno*-1);
			else scfs[workctg].push_front(vtgt[i].ctgno);
		}
	}
	int last = vtgt.size() - 1;
	if(OriVec[vtgt[last].ctgno].oristr.size() +  vtgt[last].gapdis< lib[FileNum].aveSpan){
		int sumlen = 0;
		bool bNeedMerge = 1;
		while(sumlen < lib[FileNum].aveSpan && last >= 0){
			sumlen += OriVec[vtgt[last].ctgno].oristr.size();
			sumlen += vtgt[last].gapdis;
			last--;
		}
		int startScfsw = scfs[workctg].size() - 2 * vtgt.size() - 1;
		last++;
		for(;last < vtgt.size(); last++){
			merge(scfs[workctg][startScfsw], scfdirs[workctg][startScfsw], vtgt[last].ctgno, vtgt[last].brevs, vtgt[last].gapdis);
			startScfsw += 2;
		}
	}
	for(int i=0; i<vtgt.size(); i++)bused[vtgt[i].ctgno] = 1;
	ofile<<endl;
}

void collectTargets(int workctg, deque<target> &vtgt, int rightedge, bool brevs, bool buseAll, bool need2verifyParterner){
	target t;
	int i = workctg;
	//right extension  unordered_map<int, vector<spt> >
	ofile<<"rightedge = "<<rightedge<<endl;
	if(rightedge == 469){
		rightedge = 469;
	}
	
	if(!brevs){
		ofile<<"candidates  with surpports = "<<vmpsptl[rightedge].size()<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptl[rightedge].begin(); it!=vmpsptl[rightedge].end(); it++){
			if(!buseAll && bused[it->first])continue;
			ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
			t.ctgno = it->first;
			t.brevs = 0;
			if(verifyStep1(it->second.vspt0, t, rightedge, workctg, false, need2verifyParterner)){
				t.brevs = 0;
				t.pdspt = &(it->second.vspt0);
				vtgt.push_back(t);
			};
			
			t.brevs = 1;
			ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
			if(verifyStep1(it->second.vspt1, t, rightedge, workctg, false, need2verifyParterner)){
				t.brevs = 1;
				t.pdspt = &(it->second.vspt1);
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
			if(verifyStep1(it->second.vspt0, t, rightedge, workctg, true, need2verifyParterner)){
				t.brevs = 1;
				t.pdspt = &(it->second.vspt0);
				vtgt.push_back(t);
			}
			ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
			t.brevs = 0;
			if(verifyStep1(it->second.vspt1, t, rightedge, workctg, true, need2verifyParterner)){
				t.brevs = 0;
				t.pdspt = &(it->second.vspt1);
				vtgt.push_back(t);
			}
		}
	}
	sort(vtgt.begin(), vtgt.end(), targetCmp);
}

void collectTargets2(int workctg, deque<target> &vtgt, int rightedge, bool brevs, bool buseAll, bool need2verifyParterner){
	target t;
	int i = workctg;
	//right extension  unordered_map<int, vector<spt> >
	ofile<<"rightedge = "<<rightedge<<endl;
	if(rightedge == 469){
		rightedge = 469;
	}
	
	if(!brevs){
		ofile<<"candidates  with surpports = "<<vmpsptl[rightedge].size()<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptl[rightedge].begin(); it!=vmpsptl[rightedge].end(); it++){
			if(!buseAll && bused[it->first])continue;
			ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
			t.ctgno = it->first;
			t.brevs = 0;
			if(verifyStep2(it->second.vspt0, t, rightedge, workctg, false, need2verifyParterner)){
				t.brevs = 0;
				t.pdspt = &(it->second.vspt0);
				vtgt.push_back(t);
			};
			
			t.brevs = 1;
			ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
			if(verifyStep2(it->second.vspt1, t, rightedge, workctg, false, need2verifyParterner)){
				t.brevs = 1;
				t.pdspt = &(it->second.vspt1);
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
			if(verifyStep2(it->second.vspt0, t, rightedge, workctg, true, need2verifyParterner)){
				t.brevs = 1;
				t.pdspt = &(it->second.vspt0);
				vtgt.push_back(t);
			}
			ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
			t.brevs = 0;
			if(verifyStep2(it->second.vspt1, t, rightedge, workctg, true , need2verifyParterner)){
				t.brevs = 0;
				t.pdspt = &(it->second.vspt1);
				vtgt.push_back(t);
			}
		}
	}
	sort(vtgt.begin(), vtgt.end(), targetCmp);
}

void buildGraph(unordered_map<int, set<int> > &g, deque<target> &vtgt){
	g[-1] = set<int>();
	for(int i=0; i < vtgt.size(); i++){
		g[i] = set<int>();
		g[-1].insert(i);
	}
	int vsz = vtgt.size();
	for(int i=0; i < vsz -1; i++){
		for(int j=i; j < vtgt.size(); j++){
			int bound = vtgt[i].gapdis + OriVec[vtgt[i].ctgno].oristr.size();
			if(bound < vtgt[j].gapdis){
				g[i].insert(j);
			}
		}
	}
}

void depthFirst(unordered_map<int, set<int> > &g, deque<target> &vtgt, deque<target> &path, int k, vector< deque<target> > &vdt){
	//cout<<k<<' ';
	if(g[k].empty()){ //find a path
		ofile<<"find a path : ";
		for(int i = 0; i<path.size(); i++)ofile<<path[i].ctgno<<' ';
		ofile<<endl;
		vdt.push_back(path);
		path.pop_back();
		return ;
	}else{
		for(set<int>::iterator it = g[k].begin(); it != g[k].end(); it++){
			path.push_back(vtgt[*it]);
			depthFirst(g,vtgt,path,*it,vdt);
		};
	}
	path.pop_back();
}

int zoneLength(target t){
	int res;
	int left = lib[FileNum].aveSpan - t.gapdis;
	if(left < 0)left = 0;
	int right = left - OriVec[t.ctgno].oristr.size();
	if(right < 0)right = 0;
	return (left - right);
}

double scoreOnePath(deque<target> &dt){
	double res = 0.0;
	int total = 0;
	for(int i = 0; i < dt.size(); i++)total += zoneLength(dt[i]);
	for(int i = 0; i < dt.size(); i++)res += dt[i].pdspt->size();
	res = res/total;
	return res;
}

int decision(vector< deque<target> > &vdt){
	vector<double> vsc;
	for(int i = 0; i < vdt.size(); i++){
		vsc.push_back(scoreOnePath(vdt[i]));
	}
	
	int max = 0;
	for(int i = 0; i < vsc.size(); i++){
		if(vsc[i] > vsc[max]) max = i;
	}
	cout<<"#####################"<<vsc.size()<<endl;
	sort(vsc.begin(), vsc.end(), greater<double>());
	for(int i = 0; i < vsc.size(); i++){
		cout<<vsc[i]<<' ';
	}
	cout<<endl;
	return max;
}

deque<target> selectPath(deque<target> &vtgt,int workctg, int rightedge, int dir){
	int prixlen = 0;
	int validpath = 0;
	vector<vector<int> > vv;
	unordered_map<int, set<int> > graph;
	deque<target> res;
	for(int i=0; i < (int)vtgt.size()-1; i++){
		int noOfn = (vtgt[i].gapdis - prixlen);
		prixlen += (noOfn + OriVec[vtgt[i].ctgno].oristr.size());
		if(prixlen > lib[FileNum].MaxSpan - lib[FileNum].std){
			cout<<"fuck "<<vtgt.size()<<' '<<rightedge<<endl;
			ofile<<"selectPath failed"<<endl;
			vtgt.clear();
			collectTargets2(workctg, vtgt, rightedge, scfdirs[workctg][scfdirs[workctg].size()-1], notUseAll, true);
			ofile<<"final candidates are: \n";
			for(int j=0;j<vtgt.size();j++){
				ofile<<vtgt[j].ctgno<<'\t'<<vtgt[j].partdif<<'\t'<<vtgt[j].gapdis<<'\t'<<vtgt[j].brevs<<'\t'<<OriVec[vtgt[j].ctgno].oristr.size()<<endl;
			}
			for(int i=0; i < (int)vtgt.size()-1; i++){
				int noOfn = (vtgt[i].gapdis - prixlen);
				prixlen += (noOfn + OriVec[vtgt[i].ctgno].oristr.size());
				if(prixlen > lib[FileNum].MaxSpan - lib[FileNum].std){
					cout<<"fuck again "<<vtgt.size()<<' '<<rightedge<<endl;
					vector< deque<target> > vdt;
					//cout<<"buildGraph start"<<endl;
					buildGraph(graph, vtgt);
					//cout<<"buildGraph end"<<endl;
					deque<target> path;
					//cout<<"depthFirst start"<<endl;
					depthFirst(graph, vtgt, path, -1, vdt);
					cout<<"depthFirst end"<<endl;
					int mx = decision(vdt);
					if(mx >= 0)return vdt[mx];
					return res;
				}
			}
			return vtgt;
		}
	}
	return vtgt;
}

bool scf(deque<target> &vtgt,int workctg, int rightedge, int dir){
	if(!vtgt.size())return false;
	deque<target> res = selectPath(vtgt, workctg, rightedge, dir);
	if(res.empty())return false;
	join(workctg, res, dir);
	return true;
}
