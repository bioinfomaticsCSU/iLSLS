#include<iostream>
#include<stdio.h>
#include<string.h>
#include<string>
#include<vector>
#include<algorithm>
#include<unordered_map>
#include<signal.h>
#include<sys/types.h>
#include<unistd.h>
using namespace std;

//************ debug information **************
#include<fstream>
#include<dirent.h>
ofstream ofile;
time_t ntm;
time_t prtm;
struct tm* current_time;
//*********************************************

const int NumOfThread=28;
const int Fixed = 1000;
char temp[1000000];
int subsize;
struct Pair_Str 
{
	string oristr;
	string revstr;
};
vector<Pair_Str>OriVec;
vector<Pair_Str>NextVec;

struct MateFinder 
{
	vector<string>v;
};
unordered_map<string,MateFinder>LtoR;//通过左边的mate查找右边的mate集合
unordered_map<string,MateFinder>RtoL;//通过右边的mate查找左边的mate集合
unordered_map<string,int>times;

const int MaxFileNum=10;
struct Lib
{
	char FilePe[18];
	char FileSe[18];
	int MinSpan;//min_insertsize
	int MaxSpan;//max_insertsize
	int MinOverLap;//K值
	int std;//standard deviation
	int ReadLength;//读数长度
	int errolp;
};
Lib lib[MaxFileNum];//存储文库信息

vector<bool> isused;
int ReadLength;
int FileNum;

string Reverse_Compliment(string ori)
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

bool cmp(Pair_Str p1,Pair_Str p2)
{
	return p1.oristr.size()>p2.oristr.size();
}

bool Is_Equal(string str1,string str2,int errortolerant)
{
	int i;
	int errornum=0;
	int t_size=str1.size();
	for(i=0;i<t_size;i++)
	{
		if(str1[i]!=str2[i]) errornum++;
		if(errornum>errortolerant) return false;
	}
	return true;
}

int Is_Equal2(string str1,string str2,int pos1,int pos2, int errortolerant)
{
	int i, matched(0);
	int errornum=0;
	int t_size=str1.size();
	for(i=t_size-1;i>=0;i--)
	{
		if(str1[i]!=str2[i]) errornum++;
		else matched++;
		if(errornum>errortolerant) return matched;
	}
	return matched;
}

bool Is_right_sub(string &str1,string &str2,int &overlapnum)
{
	int pos1=str1.size()-subsize;
	if(pos1<0)return false;
	string sub1=str1.substr(pos1);
	int pos2=0;
	while(1)
	{
		pos2=str2.find(sub1,pos2);
		if(pos2==-1) return false;
		if(str1.size()<pos2+subsize)
		{
			pos2+=1;
			continue;
		}
		//if(Is_Equal(str1.substr(str1.size()-pos2-subsize),str2.substr(0,pos2+subsize),1)==true) // 344801-354903   0-10102
		int matched = 0;
		matched = Is_Equal2(str1.substr(str1.size()-pos2-subsize),str2.substr(0,pos2+subsize),pos1,pos2,1);
		if(matched == pos2+subsize){
			overlapnum=pos2+subsize;
			return true;
		}else{
			if(matched > Fixed){
				str2 = str2.substr(pos2+subsize-matched);
				overlapnum=matched;
				return true;
			}
		}
		pos2+=1;
	}
	return false;
}

bool Is_left_sub(string str1,string str2,int &overlapnum)
{
	string sub1=str1.substr(0,subsize);
	int pos2=0;
	while(1)
	{
		pos2=str2.find(sub1,pos2);
		if(pos2==-1) return false;
		if(Is_Equal(str1.substr(0,str2.size()-pos2),str2.substr(pos2),1)==true)
		{
			overlapnum=str2.size()-pos2;
			return true;
		}
		pos2+=1;
	}
	return false;
}

void LeftCnt(string LeftSubStr,string RightFindStr,int & cnt)
{
	int i,j;
	unordered_map<string,MateFinder>::iterator it;
	for(i=0;i<=(int)LeftSubStr.size()-ReadLength;i++)
	{
		string sub=LeftSubStr.substr(i,ReadLength);
		it=LtoR.find(sub);
		if(it==LtoR.end()) continue;
		for(j=0;j<LtoR[sub].v.size();j++)
		{
			string tt=LtoR[sub].v[j];
			int pos=RightFindStr.find(tt);
			if(pos!=-1)
			{
				cnt++;
				break;
			}
		}
	}
	if(RightFindStr.size()<(lib[FileNum].MinSpan+lib[FileNum].MaxSpan)/2) cnt+=10;
}

void RightCnt(string RightSubStr,string LeftFindStr,int & cnt)
{
	int i,j;
	unordered_map<string,MateFinder>::iterator it;
	for(i=0;i<=(int)RightSubStr.size()-ReadLength;i++)
	{
		string sub=RightSubStr.substr(i,ReadLength);
		it=RtoL.find(sub);
		if(it==RtoL.end()) continue;
		for(j=0;j<RtoL[sub].v.size();j++)
		{
			string tt=RtoL[sub].v[j];
			int pos=LeftFindStr.find(tt);
			if(pos!=-1)
			{
				cnt++;
				break;
			}
		}
	}
	if(LeftFindStr.size()<(lib[FileNum].MinSpan+lib[FileNum].MaxSpan)/2) cnt+=10;
}

void CalculateSupport(string str1,string str2,int overlapnum,int &cnt1,int &cnt2,int &cnt3,int &cnt4,int &cnt5,int &cnt6)//str1在左，str2在右
{
	string LeftSubStr1;
	string LeftSubStr2;
	string LeftSubStr3;
	string RightFindStr;
	string RightSubStr1;
	string RightSubStr2;
	string RightSubStr3;
	string LeftFindStr;
	int i,j;
	cnt1=0; cnt2=0; cnt3=0; cnt4=0; cnt5=0; cnt6=0;
	if(str1.size()>3*ReadLength)
	{
		LeftSubStr1=str1.substr(str1.size()-3*ReadLength);
		if(str1.size()>350+3*ReadLength)
		{
			LeftSubStr2=str1.substr((str1.size()-350-3*ReadLength),3*ReadLength);
			if(str1.size()>750+3*ReadLength)
			{
				LeftSubStr3=str1.substr((str1.size()-750-3*ReadLength),3*ReadLength);
			}
			else
			{
				LeftSubStr3=str1.substr(0,3*ReadLength);
			}
		}
		else
		{
			LeftSubStr2=str1.substr(0,3*ReadLength);
			LeftSubStr3=str1.substr(0,3*ReadLength);
		}
	}
	else
	{
		LeftSubStr1=str1;
		LeftSubStr2=str1;
		LeftSubStr3=str1;
	}
	
	if(str2.size()>3*ReadLength)
	{
		RightSubStr1=str2.substr(0,3*ReadLength);
		if(str2.size()>350+3*ReadLength)
		{
			RightSubStr2=str2.substr(350,3*ReadLength);
			if(str2.size()>700+3*ReadLength)
			{
				RightSubStr3=str2.substr(700,3*ReadLength);
			}
			else
			{
				RightSubStr3=str2.substr(str2.size()-3*ReadLength);
			}
		}
		else
		{
			RightSubStr2=str2.substr(str2.size()-3*ReadLength);
			RightSubStr3=str2.substr(str2.size()-3*ReadLength);
		}
	}
	else
	{
		RightSubStr1=str2;
		RightSubStr2=str2;
		RightSubStr3=str2;
	}
	LeftFindStr=(str1.size()>lib[FileNum].MaxSpan)?str1.substr(str1.size()-lib[FileNum].MaxSpan):str1;
	RightFindStr=(str2.size()>lib[FileNum].MaxSpan)?str2.substr(0,lib[FileNum].MaxSpan):str2;
	
	LeftCnt(LeftSubStr1,RightFindStr,cnt1);
	if(cnt1==0) return;
	LeftCnt(LeftSubStr2,RightFindStr,cnt2);
	if(cnt2==0) return;
	LeftCnt(LeftSubStr3,RightFindStr,cnt3);
	if(cnt3==0) return;
	RightCnt(RightSubStr1,LeftFindStr,cnt4);
	if(cnt4==0) return;
	RightCnt(RightSubStr2,LeftFindStr,cnt5);
	if(cnt5==0) return;
	RightCnt(RightSubStr3,LeftFindStr,cnt6);
}

int main()
{
	//************** debug information ***************
	DIR *dir=NULL;
	dir = opendir("log");
	if(!dir){
		system("mkdir log");
	}

	ofile.open("log/combine.log");
	prtm = time(NULL); 
	ntm = prtm;
	current_time = localtime(&ntm);
	ofile<<"program start...\t"<< ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	//************************************************
	int i,j,k;
	int overlapnum;
	int FileId;
	int pos;
	int cnt;
	Pair_Str ps;
	Pair_Str psf;
	char SRead[110];
	if(!freopen("lib.info","r",stdin)){
		perror("lib.info");
	};
	scanf("%d",&FileNum);
	scanf("%d",&lib[0].MinOverLap);
	for(i=1;i<=FileNum;i++)
	{
		scanf("%s%d%d%d%d%d",lib[i].FilePe,&lib[i].MinSpan,&lib[i].MaxSpan,&lib[i].std,&lib[i].ReadLength, &lib[i].errolp);
		sprintf(lib[i].FilePe,"pe_%d.txt",i);
		sprintf(lib[i].FileSe,"se_%d.txt",i);
		lib[i].MinOverLap = lib[0].MinOverLap;
	}
	subsize=lib[1].MinSpan/5;

	char FileName[100];
	sprintf(FileName,"pe_%d.txt",FileNum);

	if(!freopen(FileName,"r",stdin)){
		perror(FileName);
	};
	MateFinder mf;
	char line[200];
	char *ptr;
	char ls[200],rs[200];//用C读入的字符串
	string tempstr,tempstr1,tempstr2;
	string subtempstr;
	string prefix,suffix;
	int cnt1,cnt2,cnt3,cnt4,cnt5,cnt6;
	ReadLength=0;
	while(gets(line))
	{
		if(line[0]=='>') continue;
		ptr=strtok(line, "\t"); strcpy(ls, ptr);
		ptr=strtok(NULL, "\t"); strcpy(rs, ptr);
		tempstr1=(string)ls;
		tempstr2=(string)rs;
		if(ReadLength==0) ReadLength=strlen(ls);

		LtoR[(string)ls].v.push_back((string)rs);
		
		RtoL[(string)rs].v.push_back((string)ls);
	}
	
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"read file1 over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	for(FileId=1;FileId<=FileNum;FileId++)
	{
		freopen(lib[FileId].FileSe,"r",stdin);//读入第一个文库的单读数文件
		while(scanf("%s",SRead)!=EOF)
		{
			tempstr=(string)SRead;
			if(lib[FileId].ReadLength==0) lib[FileId].ReadLength=tempstr.size();
			for(i=0;i<=lib[FileId].ReadLength-lib[1].MinOverLap;i++)
			{
				subtempstr=tempstr.substr(i,lib[1].MinOverLap);
				times[subtempstr]++;
			}
		}
	}

	if(!freopen("contigs.init.fasta","r",stdin)){
		perror("contigs.init.fasta");
	};
	
	
	while(scanf("%s",temp)!=EOF)
	{
        scanf("%s",temp);
		ps.oristr=(string)temp;
		ps.revstr=Reverse_Compliment(ps.oristr);
		OriVec.push_back(ps);
	}
	
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"read file2 over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	/*delete duplication*/
	sort(OriVec.begin(),OriVec.end(),cmp);
	isused.resize(OriVec.size());
	isused.shrink_to_fit();
	OriVec.shrink_to_fit();
	for(i=0;i<OriVec.size();i++)
	{
		if(isused[i]==true) continue;
		ps=OriVec[i];
		for(j=i+1;j<OriVec.size();j++)
		{
			if(isused[j]==true) continue;
			psf=OriVec[j];
			pos=ps.oristr.find(psf.oristr);
			if(pos!=-1)
			{
				isused[j]=true;
				break;
			}
			pos=ps.oristr.find(psf.revstr);
			if(pos!=-1)
			{
				isused[j]=true;
				break;
			}
		}
	}
	
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"delete duplication over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	for(i=0;i<OriVec.size();i++)
	{
		if(isused[i]==true) continue;
		tempstr=OriVec[i].oristr;
		prefix=tempstr.substr(0,lib[1].MinOverLap);
		if(times[prefix]<10) tempstr=tempstr.substr((tempstr.size()>lib[1].MinOverLap)?lib[1].MinOverLap:tempstr.size());
		suffix=tempstr.substr((tempstr.size()>lib[1].MinOverLap)?tempstr.size()-lib[1].MinOverLap:tempstr.size());
		if(times[suffix]<10) tempstr=tempstr.substr(0,tempstr.size()-lib[1].MinOverLap);
		OriVec[i].oristr=tempstr;
		OriVec[i].revstr=Reverse_Compliment(tempstr);
		NextVec.push_back(OriVec[i]);
	}
	
	isused.assign(isused.size(),0);
	NextVec.shrink_to_fit();
	for(i=0;i<NextVec.size();i++)
	{
		if(isused[i]==true) continue;
		while(1)
		{
			ps=NextVec[i];
			for(j=i+1;j<NextVec.size();j++)
			{
				if(isused[j]==true) continue;
				psf=NextVec[j];
				if(Is_left_sub(ps.oristr,psf.oristr,overlapnum)==true)
				{
					if(overlapnum>Fixed)
					{
						if(overlapnum>ps.oristr.size())overlapnum=ps.oristr.size();
						NextVec[i].oristr=psf.oristr+ps.oristr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
					CalculateSupport(psf.oristr,ps.oristr,overlapnum,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6);				
					if(cnt1>2&&cnt2>2&&cnt3>2&&cnt4>2&&cnt5>2&&cnt6>2)
					{
						if(overlapnum>ps.oristr.size())overlapnum=ps.oristr.size();
						NextVec[i].oristr=psf.oristr+ps.oristr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					} 
				}
				if(Is_right_sub(ps.oristr,psf.oristr,overlapnum)==true)
				{
					if(overlapnum>Fixed)
					{
						if(overlapnum>psf.oristr.size())overlapnum=psf.oristr.size();
						NextVec[i].oristr=ps.oristr+psf.oristr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
					CalculateSupport(ps.oristr,psf.oristr,overlapnum,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6);				
					if(cnt1>2&&cnt2>2&&cnt3>2&&cnt4>2&&cnt5>2&&cnt6>2)
					{
						if(overlapnum>psf.oristr.size())overlapnum=psf.oristr.size();
						NextVec[i].oristr=ps.oristr+psf.oristr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
				}
				if(Is_left_sub(ps.oristr,psf.revstr,overlapnum)==true)
				{
					if(overlapnum>Fixed)
					{
						if(overlapnum>ps.oristr.size())overlapnum=ps.oristr.size();
						NextVec[i].oristr=psf.revstr+ps.oristr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
					CalculateSupport(psf.revstr,ps.oristr,overlapnum,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6);				
					if(cnt1>2&&cnt2>2&&cnt3>2&&cnt4>2&&cnt5>2&&cnt6>2)
					{
						if(overlapnum>ps.oristr.size())overlapnum=ps.oristr.size();
						NextVec[i].oristr=psf.revstr+ps.oristr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
				}
				if(Is_right_sub(ps.oristr,psf.revstr,overlapnum)==true)
				{
					if(overlapnum>Fixed)
					{
						if(overlapnum>psf.oristr.size())overlapnum=psf.oristr.size();
						NextVec[i].oristr=ps.oristr+psf.revstr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
					CalculateSupport(ps.oristr,psf.revstr,overlapnum,cnt1,cnt2,cnt3,cnt4,cnt5,cnt6);				
					if(cnt1>2&&cnt2>2&&cnt3>2&&cnt4>2&&cnt5>2&&cnt6>2)
					{
						if(overlapnum>psf.oristr.size())overlapnum=psf.oristr.size();
						NextVec[i].oristr=ps.oristr+psf.revstr.substr(overlapnum);
						NextVec[i].revstr=Reverse_Compliment(NextVec[i].oristr);
						isused[j]=true;
						break;
					}
				}
			}
			if(j==NextVec.size()) break;
		}
	}

	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"combine over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	freopen("combine.fasta","w",stdout);
	cnt=0;
	for(i=0;i<NextVec.size();i++)
	{
		if(isused[i]==true) continue;
		cout<<">"<<cnt++<<endl<<NextVec[i].oristr.c_str()<<endl;
	}
	
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"output over,\t"<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	//************************************************
	
	_exit(0);
	return 0;
}

