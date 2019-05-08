#include<iostream>
#include<stdio.h>
#include<string.h>
#include<string>
#include<vector>
#include<stack>
#include<algorithm>
#include<cmath>
#include<unordered_map>
//************** debug information ***************
#include<fstream>
#include<dirent.h>
//************************************************
using namespace std;
char temp[1000000];
const int plusvalue=10;
const int MaxFileNum=10;
const int MaxReadLength=220;
const double NegativeMin=-1000000;
int MinExtensionLength;//k-mer扩展被接受的最小长度
int cutsize;
int FileNum;//文库的数目
int cnt=0;//k-mer的数目
struct MateFinder 
{
	vector<string>v;
	MateFinder(){};
};
unordered_map<string,MateFinder>LtoR[256*256*16];//通过左边的mate查找右边的mate集合
unordered_map<string,MateFinder>RtoL[256*256*16];//通过右边的mate查找左边的mate集合
struct Lib
{
	char FilePe[18];
	char FileSe[18];
	int MinSpan;//min_insertsize
	int MaxSpan;//max_insertsize
	int MinOverLap;//K值
	int ReadLength;//读数长度
	int std;//standard deviation
	int errolp;
};
Lib lib[MaxFileNum];//存储文库信息
struct Pair_Str 
{
	string oristr;
	string revstr;
};
vector<Pair_Str>OriVec;
vector<int>qsum;
int gapdis;//已经得到的gap距离
struct Node
{
	string CurStr;
	int index;
	int cntzero;
	double score;
	string SupportStr;
};
stack<Node>st;//保存深度优先搜索的节点
vector<Node>tempvec;//保存临时的可能节点
char line[MaxReadLength];//把输入文件按行读取然后分割
vector<bool> isdelete;
struct Element
{
	string str;
	int score;
	int maxgapsize;
	int index;
	bool breversed;
};
vector<Element>vect;
Element tempet;

//************** debug information ***************
ofstream ofile;
time_t ntm;
time_t prtm;
int indexlimit=10;
vector<unordered_map<string,int> > vumpstr2int1; //no reversed right
vector<unordered_map<string,int> > vumpstr2int2; //no reversed left
vector<unordered_map<string,int> > vumpstr2int3; //reversed right
vector<unordered_map<string,int> > vumpstr2int4; //reversed left
//************************************************

int getIndex(string str){
    int index=0;
    for(int i=0;i<str.size();i++){
        switch(str[i]){
        case 'A':
            index = index*4;
            break;
        case 'T':
            index = index*4+1;
            break;
        case 'G':
            index = index*4+2;
            break;
        case 'C':
            index = index*4+3;
            break;
        default:
            return -1;
        }
    }
    return index;
}

void findCutStr(const string& str1,const string& str2,vector<string> &vstr)
{
	string LeftSubStr1;
	string LeftSubStr2;
	string LeftSubStr3;
	string RightSubStr1;
	string RightSubStr2;
	string RightSubStr3;
	if(str1.size()>3*lib[FileNum].ReadLength)
	{
		LeftSubStr1=str1.substr(str1.size()-3*lib[FileNum].ReadLength);
		if(str1.size()>350+3*lib[FileNum].ReadLength)
		{
			LeftSubStr2=str1.substr((str1.size()-350-3*lib[FileNum].ReadLength),3*lib[FileNum].ReadLength);
			if(str1.size()>800+3*lib[FileNum].ReadLength)
			{
				LeftSubStr3=str1.substr((str1.size()-800-3*lib[FileNum].ReadLength),3*lib[FileNum].ReadLength);
			}
			else
			{
				LeftSubStr3=str1.substr(0,3*lib[FileNum].ReadLength);
			}
		}
		else
		{
			LeftSubStr2=str1.substr(0,3*lib[FileNum].ReadLength);
			LeftSubStr3=str1.substr(0,3*lib[FileNum].ReadLength);
		}
	}
	else
	{
		LeftSubStr1=str1;
		LeftSubStr2=str1;
		LeftSubStr3=str1;
	}
	if(str2.size()>3*lib[FileNum].ReadLength)
	{
		RightSubStr1=str2.substr(0,3*lib[FileNum].ReadLength);
		if(str2.size()>350+3*lib[FileNum].ReadLength)
		{
			RightSubStr2=str2.substr(350,3*lib[FileNum].ReadLength);
			if(str2.size()>800+3*lib[FileNum].ReadLength)
			{
				RightSubStr3=str2.substr(800,3*lib[FileNum].ReadLength);
			}
			else
			{
				RightSubStr3=str2.substr(str2.size()-3*lib[FileNum].ReadLength);
			}
		}
		else
		{
			RightSubStr2=str2.substr(str2.size()-3*lib[FileNum].ReadLength);
			RightSubStr3=str2.substr(str2.size()-3*lib[FileNum].ReadLength);
		}
	}
	else
	{
		RightSubStr1=str2;
		RightSubStr2=str2;
		RightSubStr3=str2;
	}
	vstr.push_back(LeftSubStr1);
	vstr.push_back(LeftSubStr2);
	vstr.push_back(LeftSubStr3);
	vstr.push_back(RightSubStr1);
	vstr.push_back(RightSubStr2);
	vstr.push_back(RightSubStr3);
}

void findLeftFindStr(string &str,string &LeftFindStr){
	LeftFindStr=(str.size()>lib[FileNum].MaxSpan)?str.substr(0,lib[FileNum].MaxSpan):str;
}

void findRightFindStr(string &str,string &RightFindStr){
	RightFindStr=(str.size()>lib[FileNum].MaxSpan)?str.substr(str.size()-lib[FileNum].MaxSpan):str;
}

void str2hash(string &str,unordered_map<string,int> & ump){
	ump.clear();
	unordered_map<string,int>::iterator it;
	for(int i=0;i<=(int)str.size()-lib[FileNum].ReadLength;i++)
	{
		string sub=str.substr(i,lib[FileNum].ReadLength);
		it = ump.find(sub);
		if(it==ump.end()){
			ump.insert(make_pair(sub,i));
		}
	}
}

bool cmp(Node n1,Node n2)
{
	if(n1.cntzero!=n2.cntzero) return n1.cntzero<n2.cntzero;
	return n1.score>n2.score;
}
bool ElmCmp(Element e1,Element e2)
{
	return e1.index<e2.index;
}

void ReadInputFile()
{
	int i;
	char *ptr;
	char ls[MaxReadLength],rs[MaxReadLength];//用C读入的字符串
	string lstr,rstr;
	int pos1,pos2;
	int index;
	unordered_map<string,MateFinder>::iterator it;
	for(i=FileNum;i<=FileNum;i++)
	{
		if(lib[i].MaxSpan<300)continue;
		if(!freopen(lib[i].FilePe,"r",stdin)){
			perror(lib[i].FilePe);
		};//读入所有文库，并保存配对关系
		while(gets(line))
		{
			if(line[0]=='>') continue;
			ptr=strtok(line, "\t"); strcpy(ls, ptr);
			ptr=strtok(NULL, "\t"); strcpy(rs, ptr);
			lstr=(string)ls;
			rstr=(string)rs;
			
			index = getIndex(lstr.substr(0,indexlimit));
			it=LtoR[index].find(lstr);
			if(it!=LtoR[index].end()) {
				(it->second).v.push_back(rstr);
			}
			index = getIndex(rstr.substr(0,indexlimit));
			it=RtoL[index].find(rstr);
			if(it!=RtoL[index].end()) {
				(it->second).v.push_back(lstr);
			}
		}
	}
	MinExtensionLength=(int)(lib[FileNum].ReadLength*1.5);
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

void LeftCnt(const string& LeftSubStr,int & cnt,int LibNum,unordered_map<string,int> & ump,int RightFindStrSize)
{
	int i,j,index;
	unordered_map<string,MateFinder>::iterator it;
	unordered_map<string,int>::iterator itump;
	for(i=0;i<=(int)LeftSubStr.size()-lib[LibNum].ReadLength;i++)
	{
		string sub=LeftSubStr.substr(i,lib[LibNum].ReadLength);
		index = getIndex(sub.substr(0,indexlimit));
		if(-1==index)continue;
		it=LtoR[index].find(sub);
		if(it==LtoR[index].end()) continue;
		MateFinder& mf=it->second;
		for(j=0;j<mf.v.size();j++)
		{
			string tt=mf.v[j];
			itump = ump.find(tt);
			if(itump!=ump.end())
			{
				if(LibNum==FileNum) 
				{
					qsum.push_back(ump[tt]+LeftSubStr.size());
				}
				cnt++;
				break;
			}
		}
	}
	if(RightFindStrSize<cutsize&&LibNum==FileNum) cnt+=plusvalue;
}

void RightCnt(const string& RightSubStr,int & cnt,int LibNum,unordered_map<string,int> & ump,int LeftFindStrSize)
{
	int i,j,index;
	unordered_map<string,MateFinder>::iterator it;
	unordered_map<string,int>::iterator itump;
	for(i=0;i<=(int)RightSubStr.size()-lib[LibNum].ReadLength;i++)
	{
		string sub=RightSubStr.substr(i,lib[LibNum].ReadLength);
		index = getIndex(sub.substr(0,indexlimit));
		if(-1==index)continue;
		it=RtoL[index].find(sub);
		if(it==RtoL[index].end()) continue;
		MateFinder& mf=it->second;
		for(j=0;j<mf.v.size();j++)
		{
			string tt=mf.v[j];
			itump = ump.find(tt);
			if(itump!=ump.end())
			{
				if(LibNum==FileNum) 
				{
					qsum.push_back(LeftFindStrSize-ump[tt]+RightSubStr.size());
				}
				cnt++;
				break;
			}
		}
	}
	if(LeftFindStrSize<cutsize&&LibNum==FileNum) cnt+=plusvalue;
}

void CalSta(const vector<int>&vec)
{
	int i;
	int sum=0;
	int numcnt=0;
	double ave=0;
	double dif=0;
	if(vec.size()==0)return;
	for(i=0;i<vec.size();i++) ave+=vec[i];
	ave=ave/vec.size();
	for(i=0;i<vec.size();i++) dif+=(vec[i]-ave)*(vec[i]-ave);
	dif=sqrt(dif/vec.size());
	for(i=0;i<vec.size();i++)
	{
		if(vec[i]>=ave-2*dif&&vec[i]<=ave+2*dif)
		{
			sum+=vec[i];
			numcnt++;
		}
	}
	if(0!=numcnt)gapdis=sum/numcnt;
}

void CalLong(const string& str1,const string& str2,int &cnt1,int &cnt2,int &cnt3,int &cnt4,int &cnt5,int &cnt6,unordered_map<string,int> & umpr,unordered_map<string,int> & umpl)
{
	string LeftSubStr1;
	string LeftSubStr2;
	string LeftSubStr3;
	string RightFindStr;
	string RightSubStr1;
	string RightSubStr2;
	string RightSubStr3;
	string LeftFindStr;
	unordered_map<string,MateFinder>::iterator it;
	cnt1=0; cnt2=0; cnt3=0; cnt4=0; cnt5=0; cnt6=0;
	if(str1.size()>3*lib[FileNum].ReadLength)
	{
		LeftSubStr1=str1.substr(str1.size()-3*lib[FileNum].ReadLength);
		if(str1.size()>350+3*lib[FileNum].ReadLength)
		{
			LeftSubStr2=str1.substr((str1.size()-350-3*lib[FileNum].ReadLength),3*lib[FileNum].ReadLength);
			if(str1.size()>800+3*lib[FileNum].ReadLength)
			{
				LeftSubStr3=str1.substr((str1.size()-800-3*lib[FileNum].ReadLength),3*lib[FileNum].ReadLength);
			}
			else
			{
				LeftSubStr3=str1.substr(0,3*lib[FileNum].ReadLength);
			}
		}
		else
		{
			LeftSubStr2=str1.substr(0,3*lib[FileNum].ReadLength);
			LeftSubStr3=str1.substr(0,3*lib[FileNum].ReadLength);
		}
	}
	else
	{
		LeftSubStr1=str1;
		LeftSubStr2=str1;
		LeftSubStr3=str1;
	}
	if(str2.size()>3*lib[FileNum].ReadLength)
	{
		RightSubStr1=str2.substr(0,3*lib[FileNum].ReadLength);
		if(str2.size()>350+3*lib[FileNum].ReadLength)
		{
			RightSubStr2=str2.substr(350,3*lib[FileNum].ReadLength);
			if(str2.size()>800+3*lib[FileNum].ReadLength)
			{
				RightSubStr3=str2.substr(800,3*lib[FileNum].ReadLength);
			}
			else
			{
				RightSubStr3=str2.substr(str2.size()-3*lib[FileNum].ReadLength);
			}
		}
		else
		{
			RightSubStr2=str2.substr(str2.size()-3*lib[FileNum].ReadLength);
			RightSubStr3=str2.substr(str2.size()-3*lib[FileNum].ReadLength);
		}
	}
	else
	{
		RightSubStr1=str2;
		RightSubStr2=str2;
		RightSubStr3=str2;
	}
	int sz1 = (str1.size()>lib[FileNum].MaxSpan)?lib[FileNum].MaxSpan:str1.size();
	int sz2 = (str2.size()>lib[FileNum].MaxSpan)?lib[FileNum].MaxSpan:str2.size();
	qsum.clear();
	LeftCnt(LeftSubStr1,cnt1,FileNum,umpl,sz2);
	RightCnt(RightSubStr1,cnt4,FileNum,umpr,sz1);
	if(cnt1<=3&&cnt4<=3)return;
	
	LeftCnt(LeftSubStr2,cnt2,FileNum,umpl,sz2);
	if((cnt1<=3||cnt2<=3)&&cnt4<=3)return;
	
	RightCnt(RightSubStr2,cnt5,FileNum,umpr,sz1);
	if((cnt1<=3||cnt2<=3)&&(cnt4<=3||cnt5<=3))return;
	
	LeftCnt(LeftSubStr3,cnt3,FileNum,umpl,sz2);
	if((cnt1<=3||cnt2<=3||cnt3<=3)&&(cnt4<=3||cnt5<=3))return;
	
	RightCnt(RightSubStr3,cnt6,FileNum,umpr,sz1);
	if((cnt1<=3||cnt2<=3||cnt3<=3)&&(cnt4<=3||cnt5<=3||cnt6<=3))return;
	
	if(qsum.size()==0) return;
	int i;
	CalSta(qsum);
}

void IsPush(const string& str,int cnt1,int cnt2,int cnt3,int cnt4,int cnt5,int cnt6,int index,bool breversed)
{
	if(str.size()<600) return;//小于600的contig之间忽略
	Element et;
	if((cnt1>3&&cnt2>3&&cnt3>3&&cnt4>3&&cnt5>3&&cnt6>3))
	{
		et.str=str;
		et.score=(cnt1)*(cnt2)*(cnt3)*(cnt4)*(cnt5)*(cnt6);
		et.maxgapsize=cutsize-gapdis;
		et.index=index;
		et.breversed=breversed;
		vect.push_back(et);
	}
	if((cnt1>1.5*plusvalue&&cnt2>1.5*plusvalue&&cnt3>1.5*plusvalue)||(cnt4>1.5*plusvalue&&cnt5>1.5*plusvalue&&cnt6>1.5*plusvalue))
	{
		if(vect.size()==0)
		{
			et.str=str;
			et.score=(cnt1)*(cnt2)*(cnt3)*(cnt4)*(cnt5)*(cnt6);
			et.maxgapsize=cutsize-gapdis;
			et.index=index;
			et.breversed=breversed;
			tempet=et;
		}
	}
}

bool IsNodeOk(const string& str)
{
	int i;
	int len=str.size();
	int scnt[4];
	memset(scnt,0,sizeof(scnt));
	for(i=0;i<len;i++)
	{
		if(str[i]=='A') scnt[0]++;
		else if(str[i]=='C') scnt[1]++;
		else if(str[i]=='G') scnt[2]++;
		else if(str[i]=='T') scnt[3]++;
	}
	for(i=0;i<4;i++) if(scnt[i]>len*0.9) return false;
	return true;
}

int GetMin(int a,int b)
{
	return a<b?a:b;
}

double GetLeftScore1(const string& AssembledString,const string& CandidateString,int index)
{
	unordered_map<string,MateFinder>::iterator it;
	int i,j;
	double score=0;
	string LeftFindStr;
	string SubStr;
	LeftFindStr=(AssembledString.size()>lib[index].MaxSpan)?AssembledString.substr(AssembledString.size()-lib[index].MaxSpan):AssembledString;
	for(i=0;i<=CandidateString.size()-lib[index].ReadLength;i++)
	{
		SubStr=CandidateString.substr(i,lib[index].ReadLength);
		index = getIndex(SubStr.substr(0,indexlimit));
		if(-1==index)continue;
		it=RtoL[index].find(SubStr);
		if(it==RtoL[index].end()) continue;
		MateFinder& mf=it->second;
		for(j=0;j<mf.v.size();j++)
		{
			string sub=mf.v[j];
			int pos=LeftFindStr.find(sub);
			if(pos==-1) continue;
			if(index==1) score+=0.1;
			else score+=0.2;
		}
	}
	return score;
}

double GetLeftScore2(const string& CandidateString,int index,const string& EndStr)
{
	unordered_map<string,MateFinder>::iterator it;
	int i,j;
	double score=0;
	string LeftFindStr;
	string SubStr;
	LeftFindStr=(EndStr.size()>lib[index].MaxSpan)?EndStr.substr(EndStr.size()-lib[index].MaxSpan):EndStr;
	for(i=0;i<=CandidateString.size()-lib[index].ReadLength;i++)
	{
		SubStr=CandidateString.substr(i,lib[index].ReadLength);
		index = getIndex(SubStr.substr(0,indexlimit));
		if(-1==index)continue;
		it=RtoL[index].find(SubStr);
		if(it==RtoL[index].end()) continue;
		MateFinder& mf=it->second;
		for(j=0;j<mf.v.size();j++)
		{
			string sub=mf.v[j];
			int pos=LeftFindStr.find(sub);
			if(pos==-1) continue;
			if(index==1) score+=0.1;
			else score+=0.2;
		}
	}
	return score;
}

void Right_Extension(string &SeedStr,string &EndStr,string &ansstr,int maxgapsize)
{
	int seedsize;
	if(maxgapsize<0) maxgapsize=0;
	seedsize=SeedStr.size();
	int nofN = maxgapsize-SeedStr.size()+seedsize;
	nofN = nofN>0?nofN:0;
	ansstr=SeedStr+string(nofN,'N')+EndStr;
}

bool PreJudge(const vector<Element>&vec)
{
	int i;
	int CalCnt=0;
	Pair_Str ps;
	for(i=0;i<vec.size();i++)
	{
		ps=OriVec[vec[i].index];
		if(ps.oristr.size()>cutsize) CalCnt++;
		if(CalCnt>1) return true;
	}
	return false;
}

bool IsAllShort(const vector<Element>&vec)
{
	int i;
	for(i=0;i<vec.size();i++)
	{
		if(vec[i].str.size()>cutsize) return false;
	}
	return true;
}

int main(int argc,char* argv[])
{
	//*************** code add by lzx *****************
	DIR *dir=NULL;
	dir = opendir("log");
	if(!dir){
		system("mkdir log");
	}
	
	ofile.open("log/scaffold.log");
	
	prtm = time(NULL); 
	ntm = prtm;
	struct tm* current_time = localtime(&ntm);
	ofile<<"program start...\t"<< ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"indexlimit = "<<indexlimit<<endl;
	//************************************************
	
	string tempstr;
	string ansstr;
	int cnt1,cnt2,cnt3,cnt4,cnt5,cnt6;
	int cntnum;
	int i,j;
	int lasttt=0;
	int times=0;
	Pair_Str ps;
	Pair_Str psf;
	if(!freopen("lib.info","r",stdin)){
		perror("lib.info");
	};//读入配置文件
	scanf("%d",&FileNum);
	scanf("%d",&lib[0].MinOverLap);
	for(i=1;i<=FileNum;i++)
	{
		scanf("%s%d%d%d%d%d",lib[i].FilePe,&lib[i].MinSpan,&lib[i].MaxSpan,&lib[i].std,&lib[i].ReadLength, &lib[i].errolp);
		sprintf(lib[i].FilePe,"pe_%d.txt",i);
		sprintf(lib[i].FileSe,"se_%d.txt",i);
	}
	lib[1].MinOverLap-=1;
	cutsize=(lib[FileNum].MinSpan+lib[FileNum].MaxSpan)/2;
	if(!freopen("combine.fasta","r",stdin)){
		perror("combine.fasta");
	};
	vector<string> vstrtmp;
	unordered_map<string,MateFinder>::iterator it;
	string LeftFindStr,RightFindStr;
	while(scanf("%s",temp)!=EOF)
	{
		vstrtmp.clear();
        scanf("%s",temp);
		ps.oristr=(string)temp;
		ps.revstr=Reverse_Compliment(ps.oristr);
		
		vumpstr2int1.push_back(unordered_map<string,int>());
		vumpstr2int2.push_back(unordered_map<string,int>());
		vumpstr2int3.push_back(unordered_map<string,int>());
		vumpstr2int4.push_back(unordered_map<string,int>());
		
		findRightFindStr(ps.oristr,RightFindStr);
		findLeftFindStr(ps.oristr,LeftFindStr);
		str2hash(RightFindStr,vumpstr2int1[OriVec.size()]);
		str2hash(LeftFindStr,vumpstr2int2[OriVec.size()]);
		
		findRightFindStr(ps.revstr,RightFindStr);
		findLeftFindStr(ps.revstr,LeftFindStr);
		str2hash(RightFindStr,vumpstr2int3[OriVec.size()]);
		str2hash(LeftFindStr,vumpstr2int4[OriVec.size()]);
		
		findCutStr(ps.oristr,ps.revstr,vstrtmp);
		findCutStr(ps.revstr,ps.oristr,vstrtmp);
		
		for(int k=0;k<vstrtmp.size();k++){
			tempstr = vstrtmp[k];
			for(i=0;i<=(int)tempstr.size()-lib[FileNum].ReadLength;i++)
			{
				string sub=tempstr.substr(i,lib[FileNum].ReadLength);
				int index = getIndex(sub.substr(0,indexlimit));
				it=LtoR[index].find(sub);
				if(it==LtoR[index].end()) {
					LtoR[index].insert(make_pair(sub,MateFinder()));
				}
				
				it=RtoL[index].find(sub);
				if(it==RtoL[index].end()) {
					RtoL[index].insert(make_pair(sub,MateFinder()));
				}
			}
		}
		OriVec.push_back(ps);
	}
	isdelete.resize(OriVec.size());
	isdelete.shrink_to_fit();
	OriVec.shrink_to_fit();
	vumpstr2int1.shrink_to_fit();
	vumpstr2int2.shrink_to_fit();
	vumpstr2int3.shrink_to_fit();
	vumpstr2int4.shrink_to_fit();
	ReadInputFile();
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"read over, "<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	prtm = ntm; 
	//************************************************
	lasttt=0;
	int nwhile = 0;
	time_t ntm2,prtm2,ntm3,prtm3;
	ntm2 = time(NULL);
	ntm3 = time(NULL);
	while(1)
	{
		nwhile++;
		for(i=0;i<OriVec.size();i++)
		{
			if(isdelete[i]==true) continue;
			cout<<"contig "<<i<<":"<<endl;
			ofile<<"contig "<<i<<":"<<endl;
			if(OriVec[i].oristr.size()<cutsize) break;
			cout<<"right extension"<<endl;
			ofile<<"right extension"<<endl;
			while(1)
			{
				ps=OriVec[i];
				vect.clear();
				tempet.maxgapsize=-1;
				for(j=i+1;j<OriVec.size();j++)
				{
					if(i==18&&j==66){
						j=66;
					}
					if(isdelete[j]==true) continue;
					psf=OriVec[j];
					CalLong(ps.oristr,psf.oristr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,vumpstr2int1[i],vumpstr2int2[j]);
					cout<<"+"<<j<<": "<<cnt1<<' '<<cnt2<<' '<<cnt3<<' '<<cnt4<<' '<<cnt5<<' '<<cnt6<<endl;
					IsPush(psf.oristr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,j,0);
					if(vect.size()>=4)break;
					CalLong(ps.oristr,psf.revstr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,vumpstr2int1[i],vumpstr2int4[j]);
					cout<<"-"<<j<<": "<<cnt1<<' '<<cnt2<<' '<<cnt3<<' '<<cnt4<<' '<<cnt5<<' '<<cnt6<<endl;
					IsPush(psf.revstr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,j,1);
					if(vect.size()>=4)break;
				}
				if(vect.size()==0&&tempet.maxgapsize!=-1) vect.push_back(tempet);
				else if(vect.size()>0)
				{
					if(IsAllShort(vect)==true&&tempet.maxgapsize!=-1) vect.push_back(tempet);
				}
				sort(vect.begin(),vect.end(),ElmCmp);
				if(vect.size()==1)
				{
					cout<<vect[0].index<<'\t'<<endl;
					ofile<<vect[0].index<<'\t'<<endl;
					Right_Extension(ps.oristr,vect[0].str,ansstr,vect[0].maxgapsize);
					OriVec[i].oristr=ansstr;
					OriVec[i].revstr=Reverse_Compliment(OriVec[i].oristr);
					isdelete[vect[0].index]=true;
				}
				else if(vect.size()>1&&vect.size()<4)
				{
					int vsize=vect.size();
					if(PreJudge(vect)==true) break;
					CalLong(vect[vsize-1].str,vect[0].str,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,\
					vect[vsize-1].breversed?vumpstr2int3[vect[vsize-1].index]:vumpstr2int1[vect[vsize-1].index],\
					vect[0].breversed?vumpstr2int4[vect[0].index]:vumpstr2int2[vect[0].index]);
					if(cnt1>3&&cnt2>3&&cnt3>3&&cnt4>3&&cnt5>3&&cnt6>3&&vect[0].str.size()>cutsize)
					{
						cout<<vect[vsize-1].index<<'\t'<<endl;
						ofile<<vect[vsize-1].index<<'\t'<<endl;
						Right_Extension(ps.oristr,vect[vsize-1].str,ansstr,vect[vsize-1].maxgapsize);
						OriVec[i].oristr=ansstr;
						OriVec[i].revstr=Reverse_Compliment(OriVec[i].oristr);
						isdelete[vect[vsize-1].index]=true;
					}
					else
					{
						cout<<vect[0].index<<'\t'<<endl;
						ofile<<vect[0].index<<'\t'<<endl;
						Right_Extension(ps.oristr,vect[0].str,ansstr,vect[0].maxgapsize);
						OriVec[i].oristr=ansstr;
						OriVec[i].revstr=Reverse_Compliment(OriVec[i].oristr);
						isdelete[vect[0].index]=true;
					}
				}
				else 
				{
					break;
				}
				findRightFindStr(OriVec[i].oristr,RightFindStr);
				str2hash(RightFindStr,vumpstr2int1[i]);
			}
			cout<<endl;
			ofile<<endl;
			cout<<"left extension"<<endl;
			ofile<<"left extension"<<endl;
			while(1)
			{
				ps=OriVec[i];
				vect.clear();
				tempet.maxgapsize=-1;
				for(j=i+1;j<OriVec.size();j++)
				{
					if(isdelete[j]==true) continue;
					psf=OriVec[j];
					CalLong(psf.oristr,ps.oristr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,vumpstr2int1[j],vumpstr2int2[i]);
					IsPush(psf.oristr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,j,0);
					if(vect.size()>=4)break;
					CalLong(psf.revstr,ps.oristr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,vumpstr2int3[j],vumpstr2int2[i]);
					IsPush(psf.revstr,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,j,1);
					if(vect.size()>=4)break;
				}
				if(vect.size()==0&&tempet.maxgapsize!=-1) vect.push_back(tempet);
				else if(vect.size()>0)
				{
					if(IsAllShort(vect)==true&&tempet.maxgapsize!=-1) vect.push_back(tempet);
				}
				sort(vect.begin(),vect.end(),ElmCmp);
				if(vect.size()==1)
				{
					cout<<vect[0].index<<'\t'<<endl;
					ofile<<vect[0].index<<'\t'<<endl;
					Right_Extension(vect[0].str,ps.oristr,ansstr,vect[0].maxgapsize);
					OriVec[i].oristr=ansstr;
					OriVec[i].revstr=Reverse_Compliment(OriVec[i].oristr);
					isdelete[vect[0].index]=true;
				}
				else if(vect.size()>1&&vect.size()<4)
				{
					int vsize=vect.size();
					if(PreJudge(vect)==true) break;
					CalLong(vect[0].str,vect[vsize-1].str,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,\
					vect[0].breversed?vumpstr2int3[vect[0].index]:vumpstr2int1[vect[0].index],\
					vect[vsize-1].breversed?vumpstr2int4[vect[vsize-1].index]:vumpstr2int2[vect[vsize-1].index]);
					if(cnt1>3&&cnt2>3&&cnt3>3&&cnt4>3&&cnt5>3&&cnt6>3&&vect[0].str.size()>cutsize)
					{
						cout<<vect[vsize-1].index<<'\t'<<endl;
						ofile<<vect[vsize-1].index<<'\t'<<endl;
						Right_Extension(vect[vsize-1].str,ps.oristr,ansstr,vect[vsize-1].maxgapsize);
						OriVec[i].oristr=ansstr;
						OriVec[i].revstr=Reverse_Compliment(OriVec[i].oristr);
						isdelete[vect[vsize-1].index]=true;
					}
					else
					{
						cout<<vect[0].index<<'\t'<<endl;
						ofile<<vect[0].index<<'\t'<<endl;
						Right_Extension(vect[0].str,ps.oristr,ansstr,vect[0].maxgapsize);
						//cout<<"-splice(i,j),  i = "<<i<<", j = "<<vect[0].index<<endl;
						OriVec[i].oristr=ansstr;
						OriVec[i].revstr=Reverse_Compliment(OriVec[i].oristr);
						isdelete[vect[0].index]=true;
					}
				}
				else 
				{
					break;
				}
				findLeftFindStr(OriVec[i].oristr,LeftFindStr);
				str2hash(LeftFindStr,vumpstr2int2[i]);
			} 
			cout<<"end"<<endl;
			//************** debug information ***************
			prtm2 = ntm2; 
			ntm2 = time(NULL);
			current_time = localtime(&ntm2);
			ofile<<"(i, nContigs) = "<<"("<<i<<", "<<OriVec.size()<<")"<<ntm2-prtm2 <<" (s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
			ofile<<"------------------------------------------"<<endl;
			//************************************************
		}
		
		int tt=0;
		for(i=0;i<OriVec.size();i++) if(isdelete[i]==true) tt++;
		//************** debug information ***************
		prtm3 = ntm3; 
		ntm3 = time(NULL);
		current_time = localtime(&ntm3);
		ofile<<"scaffold and gap round "<<nwhile<<" over,\tscaffold success "<<tt-lasttt<<" times, "<<ntm3-prtm3 <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
		ofile<<"------------------------------------------"<<endl;
		//************************************************
		
		if(lasttt!=tt) lasttt=tt;
		else break;
		for(i=0;i<OriVec.size();i++)
		{
			if(isdelete[i]==true) continue;
			findRightFindStr(OriVec[i].oristr,RightFindStr);
			findLeftFindStr(OriVec[i].oristr,LeftFindStr);
			str2hash(RightFindStr,vumpstr2int1[i]);
			str2hash(LeftFindStr,vumpstr2int2[i]);
			
			findRightFindStr(OriVec[i].revstr,RightFindStr);
			findLeftFindStr(OriVec[i].revstr,LeftFindStr);
			str2hash(RightFindStr,vumpstr2int3[i]);
			str2hash(LeftFindStr,vumpstr2int4[i]);
		}
		
		times++;
	}
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"scaffold and gap filling over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	freopen("scaffold.fasta","w",stdout);
	for(i=0;i<OriVec.size();i++)
	{
		if(isdelete[i]==true) continue;
		printf(">%d\n%s\n",i,OriVec[i].oristr.c_str());
	}
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"output over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	return 0;
}
