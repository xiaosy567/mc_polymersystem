//cmc.cpp

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#define HAVE_NO_VARIABLE_RETURN_TYPE_SUPPORT 1

//#include <mpi.h>

#include "cmc.h"


/////////////////////////////////////////////////////
void readparameter(const string TempString, const string ParaName, bool &BPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		if( TempVector[1]==string("true") || TempVector[1]==string("TRUE") ) {
			BPara=true;
        } else if( TempVector[1]==string("false") || TempVector[1]==string("FALSE") ) {
			BPara=false;
        } else {
			ErrorMSG(string(" bool variabes can only be true|TRUE|false|FALSE ... "));
			exit(LOGICERROR);
        }
	}
}
void readparameter(const string TempString, const string ParaName, int &IPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		IPara=atoi(TempVector[1].c_str());
	}
}
void readparameter(const string TempString, const string ParaName, float &FPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		FPara=atof(TempVector[1].c_str());
	}
}
void readparameter(const string TempString, const string ParaName, double &DPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		DPara=atof(TempVector[1].c_str());
	}
}
void readparameter(const string TempString, const string ParaName, string &SPara) {
	vector<string> TempVector;
	TempVector=Split(BFilter(TempString));
	if( TempVector[0] == ParaName ) {
		SPara=TempVector[1];
	}
}
///////////////////////////
cmc::cmc() {		
	ee=0.0; ff=0.0; dd=0.0;
	dd_sqrt=0.0;
	tc=NULL; ec=NULL; fc=NULL; dc=NULL;
    vecDC=NULL; coorRES=NULL;
	factor=0.0;
	tl=0.0;

	_PI_double=acos(-1.0)*2.0; _PI_half=acos(-1.0)/2.0;
	sin_PHI=0.0; cos_PHI=0.0;
	sin_THETA=0.0; cos_THETA=0.0;
	Len_CB=0.0; Len_CB_xy=0.0; Len_CB_xy_2=0.0;
	sin_Omega=0.0; cos_Omega=0.0;
	
	Omega=0.0; THETA=0.0;
	uc=0.0; vc=0.0; wc=0.0;
	TempBX=0.0; TempBY=0.0; TempBZ=0.0;
	_NUM_atoms=0; _NUM_chains=0;

	_ENER_LJ_chosenatom_backup=NULL;
	_ENER_total=0.0;
	_ENER_delta=0.0;

	_Enery_initialized=false;
	_Memory_freed=true;

	_NEIGHBOR_lj=1;

	_TEMPERATURE=1.0;
	_TEMPERATURE_REP=1.0;

	_INDEX_chn_ind=0;        // for the chosen one;
	_INDEX_chn_or_not=false; // for the chosen one;
	//_INDEX_atm_in_chn_ind=0;     // for the chosen one;
	_TYPE_atom_ind=0;        // for the chosen one;
	
	_SIZE_of_chosen_chn=0;      // for the chosen one;
	//_INDEX_chnhead=0;	_INDEX_chntail=0;

	_INDEX_chn_atom=NULL;
	_SIZE_of_chn=NULL;
	_INDEX_CHN_HEAD=NULL;	_INDEX_CHN_TAIL=NULL;
	_INDEX_chn_or_not_atom=NULL; 	_INDEX_chn_or_not_chain=NULL;
	_INDEX_atm_in_chn=NULL;   
	_TYPE_atom=NULL; 
	_INDEX_LNEIGHBOR=NULL;	_INDEX_RNEIGHBOR=NULL;
	_Index_lneighbor_ind=0;	_Index_rneighbor_ind=0;
	_INDEX_RES_ATM=NULL;       

	_BF_flag=false; _BOND_length=1.0; _BOND_delta=0.0; _PARA_K=0.0;
	_FLAG_rh2=false;
	_FLAG_rg2=false;
	stat_i=0; stat_j=0; stat_k=0; stat_size=0; stat_head=0; stat_tail=0; 
	stat_com_x=0.0; stat_com_y=0.0; stat_com_z=0.0; tempnum=0.0;
	_SIZE_memo=0;
	_IFVERBOSE=false;

	_RUNTIMES_eachstep=0;         // 100000
    _RUNTIMES_totalnum=0;         // 10000
	_RUNTIMES_output=0;           // 40 

	_I_eachstep=0; _I_totalnum=0;
	
	iy=0; jj=0; ran=0.0; ir=NULL;
	iseed_zero=-1; iseed_len1=-1; iseed_angle=-1; iseed_index=-1; iseed_rand=-1;

	//_RUNTIMES_iteration=0;
	//_START_fromzero=false;
	//_ITER=0;

	ddelta=0.0; rmin=0.0; rmax=0.0; r2min=0.0; r2max=0.0;
	
	//_DIS=0.0;
	_DIS2=0.0;_DIS3=0.0;_DIS6=0.0;_DIS9=0.0;_DIS12=0.0;
	_DIS_x=0.0; _DIS_y=0.0; _DIS_z=0.0;
	tempind=0; tempind_x=0; tempind_y=0;

	//tempindex_judge=0; tempindex_judge_bak=0;
	tempnumer_ret=0; TempEner=0.0;

	_XX=NULL; _YY=NULL; _ZZ=NULL;
	//_XX_rec=NULL; _YY_rec=NULL; _ZZ_rec=NULL;

	_MC_NUM_TOT=0;
	_MC_NUM_SUC=0;
	_MC_NUM_FIL=0;
	_Statistic_over=false;
	//_MC_NUM_STT=0;
}
///////////////////////////
cmc::~cmc() {
	this->memo_free();
}
///////////////////////////
void cmc::memo_allocation() {	
	//int i=0;
	string TempStr;
	string TempInterval;

	TempStr=string(" [ memory allocating... ");
	TempInterval=string("-");
	cout<<TempStr;
	///////////////   code start  //////////////////	
	ec=new double[3];
	fc=new double[3];
	dc=new double[3];
	tc=new double[3];

	vecDC=new double[3];
	coorRES=new double[3];
	
	ir=new int[98];
	
	if(_NUM_atoms==0) {
		ErrorMSG(string("_NUM_atoms==0;::cmc::memo_allocation"));
		exit(SIZEERROR);
	}
	_SIZE_memo=MEMOSIZE(_NUM_atoms);
	_ENER_LJ_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_ENER_LJ_chosenatom_backup=new double[_SIZE_memo];

	_ENER_BF_eachatom=new double[_SIZE_memo];
	_ENER_BF_chosenatom_backup=0.0;

	_DIS2_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS2_chosenatom_backup=new double[_SIZE_memo];
	_DIS_x_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS_x_chosenatom_backup=new double[_SIZE_memo];
	_DIS_y_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS_y_chosenatom_backup=new double[_SIZE_memo];
	_DIS_z_eachatom.Build(_SIZE_memo,_SIZE_memo);
	_DIS_z_chosenatom_backup=new double[_SIZE_memo];
	_COM_x=new double[_NUM_chains];
	_COM_y=new double[_NUM_chains];
	_COM_z=new double[_NUM_chains];
	_RG2_x=new double[_NUM_chains];
	_RG2_y=new double[_NUM_chains];
	_RG2_z=new double[_NUM_chains];
	_rgyration2=new double[_NUM_chains];
	_rhead2=new double[_NUM_chains];

	_INDEX_chn_atom=new int[_SIZE_memo];
	_SIZE_of_chn=new int[_NUM_chains];
	_INDEX_CHN_HEAD=new int[_NUM_chains];
	_INDEX_CHN_TAIL=new int[_NUM_chains];
	_INDEX_atm_in_chn=new int[_SIZE_memo];  
	_TYPE_atom=new int[_SIZE_memo];
	_INDEX_LNEIGHBOR=new int[_SIZE_memo];
	_INDEX_RNEIGHBOR=new int[_SIZE_memo];
	_INDEX_RES_ATM=new int[_SIZE_memo];
    _XX=new double[_SIZE_memo];
	_YY=new double[_SIZE_memo];
	_ZZ=new double[_SIZE_memo];
	/*_NUM_rec=_SIZE_memo*_RUNTIMES_eachstep;
	_XX_rec=new double[_NUM_rec];
	_YY_rec=new double[_NUM_rec];
	_ZZ_rec=new double[_NUM_rec];*/

	_INDEX_chn_or_not_atom=new int[_SIZE_memo];
	_INDEX_chn_or_not_chain=new int[_NUM_chains];

	_EPSILON_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_SIGMA_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_SIGMA3_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_SIGMA6_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_SIGMA9_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_SIGMA12_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_E_cut_RR_eachatom.Build(_SIZE_memo, _SIZE_memo);
	_R_cut_RR_eachatom.Build(_SIZE_memo, _SIZE_memo);

    ///////////////   code start  //////////////////
    cout<<setw(_OUTPUT_LEN-TempStr.size())<<TempInterval<<" done! ]"<<endl;
	_Memory_freed=false;
}
////////////////////////////
void cmc::memo_setzero() {
	
	string TempStr;
	string TempInterval;
	
	TempStr=string(" [ memory resetting ");
	TempInterval=string("-");
	cout<<TempStr;
	///////////////   code start  //////////////////
	MEMOSETZERO(ec, sizeof(double)*3);
	MEMOSETZERO(fc, sizeof(double)*3);
	MEMOSETZERO(dc, sizeof(double)*3);

	MEMOSETZERO(coorRES, sizeof(double)*3);
	MEMOSETZERO(vecDC, sizeof(double)*3);

	MEMOSETZERO(tc, sizeof(double)*3);

	MEMOSETZERO(ir, sizeof(int)*98);

	_ENER_LJ_eachatom.SetZero();
	MEMOSETZERO(_ENER_LJ_chosenatom_backup, sizeof(double)*_SIZE_memo);

	MEMOSETZERO(_ENER_BF_eachatom, sizeof(double)*_SIZE_memo);

	_DIS2_eachatom.SetZero();
	MEMOSETZERO(_DIS2_chosenatom_backup, sizeof(double)*_SIZE_memo);
	_DIS_x_eachatom.SetZero();
	MEMOSETZERO(_DIS_x_chosenatom_backup, sizeof(double)*_SIZE_memo);
	_DIS_y_eachatom.SetZero();
	MEMOSETZERO(_DIS_y_chosenatom_backup, sizeof(double)*_SIZE_memo);
	_DIS_z_eachatom.SetZero();
	MEMOSETZERO(_DIS_z_chosenatom_backup, sizeof(double)*_SIZE_memo);

	MEMOSETZERO(_COM_x, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_COM_y, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_COM_z, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_x, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_y, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_RG2_z, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_rgyration2, sizeof(double)*_NUM_chains);
	MEMOSETZERO(_rhead2, sizeof(double)*_NUM_chains);

	MEMOSETZERO(_INDEX_chn_atom, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_SIZE_of_chn, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_CHN_HEAD, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_CHN_TAIL, sizeof(int)*_NUM_chains);
	MEMOSETZERO(_INDEX_atm_in_chn, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_TYPE_atom, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_LNEIGHBOR, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_RNEIGHBOR, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_RES_ATM, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_XX, sizeof(double)*_SIZE_memo);
	MEMOSETZERO(_YY, sizeof(double)*_SIZE_memo);
	MEMOSETZERO(_ZZ, sizeof(double)*_SIZE_memo);
	/*MEMOSETZERO(_XX_rec, sizeof(double)*_NUM_rec);
	MEMOSETZERO(_YY_rec, sizeof(double)*_NUM_rec);
	MEMOSETZERO(_ZZ_rec, sizeof(double)*_NUM_rec);*/

	MEMOSETZERO(_INDEX_chn_or_not_atom, sizeof(int)*_SIZE_memo);
	MEMOSETZERO(_INDEX_chn_or_not_chain, sizeof(int)*_NUM_chains);

	_EPSILON_eachatom.SetZero();
	_SIGMA_eachatom.SetZero();
	_SIGMA3_eachatom.SetZero();
	_SIGMA6_eachatom.SetZero();
	_SIGMA9_eachatom.SetZero();
	_SIGMA12_eachatom.SetZero();
	_E_cut_RR_eachatom.SetZero();
	_R_cut_RR_eachatom.SetZero();

    ///////////////   code start  //////////////////
    cout<<setw(_OUTPUT_LEN-TempStr.size())<<TempInterval<<" done! ]"<<endl;
}
///////////////////////////
void cmc::memo_evaluation() {
	//string TempStr_nouse;
	string TempStr;
	string TempInterval;
	
	TempStr=string(" [ memory evaluating ");
	TempInterval=string("-");
	cout<<TempStr;

	int Size_RES=0;
	int Size_ATM=0;
	int index_atm=0;
	int index_res=0;

	int i=0;
	int j=0;
	int k=0;
	//int l=0;
	int numerator=0;
	/*for(i=1; i<=_NUM_atoms; i++) {
		_INDEX_chn_atom[i]=_system_._Index_CHN[i];
		_INDEX_atm_in_chn[i]=_system_._Index_ATM_in_CHN[i];
		_TYPE_atom[i]=_system_._Type_ATOM[i];
		_XX[i]=_system_._XX[i];
		_YY[i]=_system_._YY[i];
		_ZZ[i]=_system_._ZZ[i];
	}*/
	int _CHCK_UPPER_POINT=20;

	for(i=0; i<_NUM_chains; i++) {
		//cout<<i<<endl;
		Size_RES=_system_.chains[i].nresidues;
		if( _system_.chains[i].index_chn_real<_CHCK_UPPER_POINT ) {
			_INDEX_chn_or_not_chain[i]=true;
		} else {
			_INDEX_chn_or_not_chain[i]=false;
		}
		index_atm=0;
		for(j=0; j<Size_RES; j++) {
			//cout<<j<<endl;
			Size_ATM=_system_.chains[i].residues[j].natoms;
			//cout<<Size_ATM<<" "<<Size_RES<<endl;
			for(k=0; k<Size_ATM; k++) {
				//cout<<"k="<<k<<endl;
				_XX[numerator+1]=_system_.chains[i].residues[j].atoms[k].x_coordinate;
				_YY[numerator+1]=_system_.chains[i].residues[j].atoms[k].y_coordinate;
				_ZZ[numerator+1]=_system_.chains[i].residues[j].atoms[k].z_coordinate;
				//cout<<_XX[numerator+1]<<" "<<_YY[numerator+1]<<" "<<_ZZ[numerator+1]<<endl; 
				if( j==0 && k==0 ) {
					_TYPE_atom[numerator+1]=-1; //head
					_INDEX_CHN_HEAD[i]=numerator+1;
				} else if( j==(Size_RES-1) && k==(Size_ATM-1) ) {
					_TYPE_atom[numerator+1]=1; //tail
					_INDEX_CHN_TAIL[i]=numerator+1;
				} else {
					_TYPE_atom[numerator+1]=0;
				}
				_INDEX_chn_atom[numerator+1]=i;
				_INDEX_chn_or_not_atom[numerator+1]=_INDEX_chn_or_not_chain[i];
				_INDEX_RES_ATM[numerator+1]=index_res+j;
				//_Index_ATM_in_RES[numerator+1]=k;
				_INDEX_atm_in_chn[numerator+1]=index_atm+k;
				numerator++;
			    //cout<<_Index_CHN_real[numerator+1]<<endl;
			}
			index_atm+=Size_ATM;
		}
		_SIZE_of_chn[i]=index_atm;
		index_res+=Size_RES;
	}
	for(i=1;i<_SIZE_memo;i++) {
		_INDEX_LNEIGHBOR[i]=(i-_NEIGHBOR_lj)<=_INDEX_CHN_HEAD[_INDEX_chn_atom[i]]?_INDEX_CHN_HEAD[_INDEX_chn_atom[i]]:(i-_NEIGHBOR_lj); //[
		_INDEX_RNEIGHBOR[i]=(i+_NEIGHBOR_lj)>=_INDEX_CHN_TAIL[_INDEX_chn_atom[i]]?_INDEX_CHN_TAIL[_INDEX_chn_atom[i]]:(i+_NEIGHBOR_lj); //]
	}
	if( _NUM_atoms==1 && _INDEX_chn_or_not_atom[1] ) {
		cout<<" _NUM_atoms="<<_NUM_atoms<<endl;
		cout<<" _INDEX_chn_or_not_atom[1]="<<_INDEX_chn_or_not_atom[1]<<endl;
		cout<<" chain with only one atom can not move ... "<<endl;
		cout<<" try again and make sure if this atom belongs to a chain or not (idx_chn > 20 ?)"<<endl;
		//exit(-1);
		exit(LOGICERROR);
	}
	//_sigma_surf_1_6=pow(0.4, 1.0/6.0)*_sigma_surf; //new_edition for 9, 3 LJ of semi-finite substrate;
	cout<<setw(_OUTPUT_LEN-TempStr.size())<<TempInterval<<" done! ]"<<endl;
	check_memo();
	///////////////   code start  //////////////////
}
//////////////////////////
inline void cmc::check_memo() {
	cout<<endl<<" check the molecule information .... "<<endl;
	int i;
	for(i=0; i<_NUM_chains; i++) {
		if(_INDEX_chn_or_not_chain[i]) {
			cout<<" the size of chain ["<<i<<"] is: "<<_SIZE_of_chn[i]<<"("<<_INDEX_CHN_HEAD[i]<<","<<_INDEX_CHN_TAIL[i]<<")"<<endl;
		} else {
			cout<<" the size of non-chain ["<<i<<"] is: "<<_SIZE_of_chn[i]<<"("<<_INDEX_CHN_HEAD[i]<<","<<_INDEX_CHN_TAIL[i]<<")"<<endl;
		}
	}
	int tnum_chn=int( log10(double(_NUM_chains)))+1;
	int tidx_atm=0;
	for(i=1; i<_SIZE_memo; i++) {
		tidx_atm=int( log10(double(_SIZE_of_chn[_INDEX_chn_atom[i]]) ))+1;
		if(_INDEX_chn_or_not_atom[i]) {
			cout<<" ccatom [";
		} else {
			cout<<" ncatom [";
		}
		cout<<setw(tidx_atm)<<i<<"] on chn:["<<setw(tnum_chn)<<_INDEX_chn_atom[i]<<"] the ["
		    <<setw(tidx_atm)<<_INDEX_atm_in_chn[i]+1<<"]th one ["
		    <<setw(tidx_atm)<<_INDEX_LNEIGHBOR[i]<<","<<setw(tidx_atm)<<_INDEX_RNEIGHBOR[i]<<"], type ["
		    <<setw(2)<<_TYPE_atom[i]
		    <<"] idx_r_in_m: ["<<_INDEX_RES_ATM[i]
		    <<"] c:["<<_XX[i]<<","<<_YY[i]<<","<<_ZZ[i]<<"];"<<endl;
	}
	cout<<" check the molecule information done! "<<endl<<endl;
}
///////////////////////////
void cmc::memo_free() {	
	cout<<" memory deallocation ... "<<endl;
	if(_Memory_freed) { return; }
	delete[] ec;
	delete[] fc;
	delete[] dc;
	delete[] coorRES;
	delete[] vecDC;
	delete[] tc;
	delete[] ir;

	_ENER_LJ_eachatom.Release();
	delete[] _ENER_LJ_chosenatom_backup;

	delete[] _ENER_BF_eachatom;

	_DIS2_eachatom.Release();
	delete[] _DIS2_chosenatom_backup;
	_DIS_x_eachatom.Release();
	delete[] _DIS_x_chosenatom_backup;
	_DIS_y_eachatom.Release();
	delete[] _DIS_y_chosenatom_backup;
	_DIS_z_eachatom.Release();
	delete[] _DIS_z_chosenatom_backup;

	delete[] _COM_x;
	delete[] _COM_y;
	delete[] _COM_z;
	delete[] _RG2_x;
	delete[] _RG2_y;
	delete[] _RG2_z;
	delete[] _rgyration2;
	delete[] _rhead2;

	delete[] _INDEX_chn_atom;
	delete[] _SIZE_of_chn;
	delete[] _INDEX_CHN_HEAD;
	delete[] _INDEX_CHN_TAIL;
	delete[] _INDEX_atm_in_chn;  
	delete[] _TYPE_atom;
	delete[] _INDEX_LNEIGHBOR;
	delete[] _INDEX_RNEIGHBOR;
	delete[] _INDEX_RES_ATM;
    delete[] _XX;
	delete[] _YY;
	delete[] _ZZ;
	/*delete[] _XX_rec;
	delete[] _YY_rec;
	delete[] _ZZ_rec;*/

	delete[] _INDEX_chn_or_not_atom;
	delete[] _INDEX_chn_or_not_chain;
	
	_EPSILON_eachatom.Release();
	_SIGMA_eachatom.Release();
	_SIGMA3_eachatom.Release();
	_SIGMA6_eachatom.Release();
	_SIGMA9_eachatom.Release();
	_SIGMA12_eachatom.Release();
	_E_cut_RR_eachatom.Release();
	_R_cut_RR_eachatom.Release();

	if(!_Statistic_over) { close_statistic(); }
	cout<<" memory deallocation done!"<<endl;
	_Memory_freed=true;
}
///////////////////////////
void cmc::load_parameters(const string FILENAME_para) {
	ifstream para_ifstream(FILENAME_para.c_str());
	if(para_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<FILENAME_para<<endl;
		exit(IOERROR);
	}
	cout<<endl<<" Loading parameters from file [ "<<FILENAME_para<<" ] ..."<<endl;
	string tempstr;
	vector<string> tempvec;
	while( getline(para_ifstream, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		if(tempvec.size()!=0) {
			readparameter(tempstr, string("_IFVERBOSE"), _IFVERBOSE);

			readparameter(tempstr, string("_FILENAME_conf"), _FILENAME_conf);
			
			readparameter(tempstr, string("_CoorX1"), _CoorX1);
			readparameter(tempstr, string("_CoorX2"), _CoorX2);
			readparameter(tempstr, string("_CoorY1"), _CoorY1);
			readparameter(tempstr, string("_CoorY2"), _CoorY2);
			readparameter(tempstr, string("_CoorZ1"), _CoorZ1);
			readparameter(tempstr, string("_CoorZ2"), _CoorZ2);
			
			readparameter(tempstr, string("_PBC_Dim"), _PBC_Dim);
			readparameter(tempstr, string("_BF_flag"), _BF_flag);
			readparameter(tempstr, string("_PARA_K"), _PARA_K);
			readparameter(tempstr, string("_BOND_delta"), _BOND_delta);
			readparameter(tempstr, string("_BOND_length"), _BOND_length);
			readparameter(tempstr, string("_TEMPERATURE"), _TEMPERATURE);
			readparameter(tempstr, string("_EPSILON_FN"), _EPSILON_FN);
			readparameter(tempstr, string("_SIGMA_FN"), _SIGMA_FN);
			readparameter(tempstr, string("_NEIGHBOR_lj"), _NEIGHBOR_lj);
			readparameter(tempstr, string("_FLAG_rh2"), _FLAG_rh2);
			readparameter(tempstr, string("_FLAG_rg2"), _FLAG_rg2);

			readparameter(tempstr, string("_RUNTIMES_eachstep"), _RUNTIMES_eachstep);
			readparameter(tempstr, string("_RUNTIMES_totalnum"), _RUNTIMES_totalnum);
			readparameter(tempstr, string("_RUNTIMES_output"), _RUNTIMES_output);

		}
	}
	if(_RUNTIMES_output>_RUNTIMES_totalnum) {
		ErrorMSG("_RUNTIMES_output>_RUNTIMES_totalnum::load_parameters");
		exit(LOGICERROR);	
	}
	if(_BOND_length<=0.0) {
		ErrorMSG("_BOND_length<=0.0::load_parameters");
		exit(LOGICERROR);		
	}
	if(_NEIGHBOR_lj<=0) {
		ErrorMSG("_NEIGHBOR_lj<=0.0::load_parameters");
		exit(LOGICERROR);
	}
	if(_CoorX1>=_CoorX2 || _CoorY1>=_CoorY2 || _CoorZ1>=_CoorZ2) {
		ErrorMSG("_CoorX1>=_CoorX2 or _CoorY1>=_CoorY2 or _CoorZ1>=_CoorZ2::load_parameters");
		exit(LOGICERROR);
	}
	if(_PBC_Dim<1) {
		ErrorMSG("_PBC_Dim<1::load_parameters");
		exit(LOGICERROR);
	}
	
	cout<<setiosflags(ios::left);

	cout<<setw(30)<<" ::_IFVERBOSE: "<<_IFVERBOSE<<endl;
	cout<<setw(30)<<" ::_FILENAME_conf: "<<_FILENAME_conf.c_str()<<endl;
	cout<<setw(30)<<" ::_CoorX1: "<<_CoorX1<<endl;
	cout<<setw(30)<<" ::_CoorX2: "<<_CoorX2<<endl;
	cout<<setw(30)<<" ::_CoorY1: "<<_CoorY1<<endl;
	cout<<setw(30)<<" ::_CoorY2: "<<_CoorY2<<endl;
	cout<<setw(30)<<" ::_CoorZ1: "<<_CoorZ1<<endl;
	cout<<setw(30)<<" ::_CoorZ2: "<<_CoorZ2<<endl;
	cout<<setw(30)<<" ::_PBC_Dim: "<<_PBC_Dim<<endl;
	_PBL_X=_CoorX2-_CoorX1;
	_PBL_Y=_CoorY2-_CoorY1;
	_PBL_Z=_CoorZ2-_CoorZ1;
	cout<<setw(30)<<" ::_PBL_X: "<<_PBL_X<<endl;	
	cout<<setw(30)<<" ::_PBL_Y: "<<_PBL_Y<<endl;
	cout<<setw(30)<<" ::_PBL_Z: "<<_PBL_Z<<endl;
	if(_BOND_length<=0.0) {
		ErrorMSG("_BOND_length<=0.0::load_parameters");
		exit(LOGICERROR);		
	}

	cout<<setw(30)<<" ::_BF_flag: "<<_BF_flag<<endl;
	cout<<setw(30)<<" ::_PARA_K: "<<_PARA_K<<endl;
	cout<<setw(30)<<" ::_BOND_delta: "<<_BOND_delta<<endl;
	cout<<setw(30)<<" ::_BOND_length: "<<_BOND_length<<endl;
	if(_BF_flag) {	
		rmin=_BOND_length-_BOND_delta;
		rmax=_BOND_length+_BOND_delta;
		ddelta=2.0*_BOND_delta;
	}
	cout<<setw(30)<<" ::rmin: "<<rmin<<endl;
	cout<<setw(30)<<" ::rmax: "<<rmax<<endl;
	cout<<setw(30)<<" ::ddelta: "<<ddelta<<endl;
	if(_TEMPERATURE<=0.0) {
		cout<<"_TEMPERATURE="<<_TEMPERATURE<<"<=0.0"<<endl;
		cout<<"exit!"<<endl;
		exit(LOGICERROR);
	} else {
		_TEMPERATURE_REP=1/_TEMPERATURE;
	}
	cout<<setw(30)<<" ::_TEMPERATURE: "<<_TEMPERATURE<<endl;
	cout<<setw(30)<<" ::_TEMPERATURE_REP: "<<_TEMPERATURE_REP<<endl;

	//make_dir();

	cout<<setw(30)<<" ::_EPSILON_FN: "<<_EPSILON_FN<<endl;
	cout<<setw(30)<<" ::_SIGMA_FN: "<<_SIGMA_FN<<endl;
	cout<<setw(30)<<" ::_NEIGHBOR_lj: "<<_NEIGHBOR_lj<<endl;
	cout<<setw(30)<<" ::_FLAG_rh2: "<<_FLAG_rh2<<endl;
	cout<<setw(30)<<" ::_FLAG_rg2: "<<_FLAG_rg2<<endl;

	//cout<<setw(30)<<" ::_RUNTIMES_eachstep: "<<_RUNTIMES_eachstep<<endl;
	cout<<setw(30)<<" ::_RUNTIMES_totalnum: "<<_RUNTIMES_totalnum<<endl;
	cout<<setw(30)<<" ::_RUNTIMES_output: "<<_RUNTIMES_output<<endl;

	//cout<<setw(30)<<" ::_RUNTIMES_iteration: "<<_RUNTIMES_iteration<<endl;
	//cout<<setw(30)<<" ::_START_fromzero: "<<_START_fromzero<<endl;
	cout<<resetiosflags(ios::left);

	para_ifstream.close();
	cout<<" Parameters loaded!"<<endl;
}
void cmc::write_parameters(const string FILENAME_para) {
	string OFName=FILENAME_para;
	int OFNameLen=FILENAME_para.size();
	OFName=OFName.substr(0, OFNameLen-4)+string(".cmp")+OFName.substr(OFNameLen-4, OFNameLen);
	ofstream para_ofstream( OFName.c_str() );
	if(para_ofstream==NULL) {
		cout<<" Error: Can not open file: "<<OFName<<endl;
		exit(IOERROR);
	}
	cout<<setiosflags(ios::left);

	writeparameter(para_ofstream, string("_IFVERBOSE"), _IFVERBOSE);

	writeparameter(para_ofstream, string("_FILENAME_conf"), _FILENAME_conf.c_str());
	writeparameter(para_ofstream, string("_CoorX1"), _CoorX1);
	writeparameter(para_ofstream, string("_CoorX2"), _CoorX2);
	writeparameter(para_ofstream, string("_CoorY1"), _CoorY1);
	writeparameter(para_ofstream, string("_CoorY2"), _CoorY2);
	writeparameter(para_ofstream, string("_CoorZ1"), _CoorZ1);
	writeparameter(para_ofstream, string("_CoorZ2"), _CoorZ2);
	writeparameter(para_ofstream, string("_PBC_Dim"), _PBC_Dim);
	writeparameter(para_ofstream, string("_PBL_X"), _PBL_X);
	writeparameter(para_ofstream, string("_PBL_Y"), _PBL_Y);
	writeparameter(para_ofstream, string("_PBL_Z"), _PBL_Z);
	writeparameter(para_ofstream, string("_BF_flag"), _BF_flag);
	writeparameter(para_ofstream, string("_PARA_K"), _PARA_K);
	writeparameter(para_ofstream, string("_BOND_delta"), _BOND_delta);
	writeparameter(para_ofstream, string("_BOND_length"), _BOND_length);
	writeparameter(para_ofstream, string("rmin"), rmin);
	writeparameter(para_ofstream, string("rmax"), rmax);
	writeparameter(para_ofstream, string("ddelta"), ddelta);

	cout<<setw(30)<<" ::_NUM_atoms: "<<_NUM_atoms<<endl;
    writeparameter(para_ofstream, string("_NUM_atoms"), _NUM_atoms);
    cout<<setw(30)<<" ::_NUM_chains: "<<_NUM_chains<<endl;
    writeparameter(para_ofstream, string("_NUM_chains"), _NUM_chains);
	writeparameter(para_ofstream, string("_TEMPERATURE_REP"), _TEMPERATURE_REP);
	writeparameter(para_ofstream, string("_TEMPERATURE"), _TEMPERATURE);
	writeparameter(para_ofstream, string("_TEMPERATURE_REP"), _TEMPERATURE_REP);

	writeparameter(para_ofstream, string("_SIGMA_FN"), _SIGMA_FN);
	writeparameter(para_ofstream, string("_NEIGHBOR_lj"), _NEIGHBOR_lj);
	writeparameter(para_ofstream, string("_FLAG_rh2"), _FLAG_rh2);
	writeparameter(para_ofstream, string("_FLAG_rg2"), _FLAG_rg2);

	//writeparameter(para_ofstream, string("_NUM_replica"), _NUM_replica);
	_RUNTIMES_eachstep*=_NUM_atoms;
	cout<<setw(30)<<" ::_RUNTIMES_eachstep: "<<_RUNTIMES_eachstep<<endl;
	writeparameter(para_ofstream, string("_RUNTIMES_eachstep"), _RUNTIMES_eachstep);
	writeparameter(para_ofstream, string("_RUNTIMES_totalnum"), _RUNTIMES_totalnum);
	writeparameter(para_ofstream, string("_RUNTIMES_output"), _RUNTIMES_output);

	//writeparameter(para_ofstream, string("_RUNTIMES_iteration"), _RUNTIMES_iteration);
	//writeparameter(para_ofstream, string("_START_fromzero"), _START_fromzero);

	cout<<resetiosflags(ios::left);
	para_ofstream.close();
	cout<<" check [ "<<OFName<<" ] to c if paras loaded right or not."<<endl;
}
///////////////////////////
/*void cmc::make_dir() {
	//using sstream::stringstream;
	//stringstream ss;
	//ss<<_TEMPERATURE;
	//string temp_dir=ss.str();
	//int temp_len=temp_dir.length();
	//_OUTPUTDIR=new char[10];
	//temp_dir.copy(_OUTPUTDIR,temp_len,0);
	//_OUTPUTDIR[temp_len]='\0';
	//sprintf(_OUTPUTDIR, "%010.5f",_TEMPERATURE);
	_FNoutput=new char[15];
	sprintf(_FNoutput, "%010.5f.dat",_TEMPERATURE);
	ifstream check_file(_FNoutput);
	if( check_file ) {
		system( (string("rm -fr ")+string(_FNoutput)).c_str() );
		//cout<<" DIR: "<<temp_dir<<" was created!"<<endl;
	} 
	check_file.close();
	//system( (string("mkdir ")+string(_OUTPUTDIR)).c_str() );
	cout<<" :: FILE: ["<<_FNoutput<<"] was created!"<<endl;
	///okay, no problem.
	//FNAME=new char[10int(_RUNTIMES_totalnum/10)+];
	//sprintf(FNAME, "%s/%020.0f.dat",_OUTPUTDIR,_MC_NUM_TOT);
}*/
///////////////////////////
void cmc::init_conformation_a(const int NUM_atoms, const int NUM_chains, const double LEN_bond, const string FILENAME_conf) {
	int i=0;
	int j=0;
	string AType=string("ATOM");
	string AChem_Name=string("C");
	string ASpec_Name=string("");
	string AResi_Name=string("CHN");
	string AChn_Name=string("A");
	char AChn_Name_temp[2];//for chain index use!
	int AResi_Index=1;
	double sqrt_3_d_2_BOND=sqrt(3)/2*LEN_bond;
	double sqrt_1_d_2_BOND=0.5*LEN_bond;
	int temp_atoms=NUM_atoms;
	for(i=0; i<NUM_chains; i++) {
		AChn_Name_temp[0]=AChn_Name.c_str()[0]+i;
		AChn_Name_temp[1]='\0';
		AChn_Name=string(AChn_Name_temp);
		for(j=0; j<temp_atoms; j++) {
			_system_.add_arbitrary_info(AType, j+i*temp_atoms+1, AChem_Name, ASpec_Name, AResi_Name, AChn_Name, AResi_Index,
										sqrt_3_d_2_BOND*j, 
										sqrt_1_d_2_BOND*(j%2), 
										i*LEN_bond*3.0+LEN_bond);
		}
	}
	_system_.writepdbinfo(FILENAME_conf.c_str(), true);
	_system_.memo_free();
	cout<<endl<<" "<<NUM_atoms*NUM_chains<<" atoms was constructed on "<<NUM_chains<<" chains!"<<endl;
}
///////////////////////////
void cmc::load_conformation() {//initialization conformation; init.pdb (xyz);
	cout<<" loading conformation from: [ "<<_FILENAME_conf<<" ];"<<endl;
	_system_.readpdbinfo(_FILENAME_conf.c_str(), 1.0, _IFVERBOSE);
	_NUM_chains=_system_.nchains;
	_NUM_atoms=0;
	int j;
	int tnresidues;
	for(int i=0; i<_NUM_chains; i++) {
		tnresidues=_system_.chains[i].nresidues;
		for(j=0; j<tnresidues; j++) {
			_NUM_atoms+=_system_.chains[i].residues[j].natoms;
		}
	}
}
///////////////////////////
void cmc::init_epsilonsigma() {//temporary arbitrary!
	ifstream epsilon_ifstream(_EPSILON_FN.c_str());
	if(epsilon_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<_EPSILON_FN<<endl;
		exit(IOERROR);
	}
	
	cout<<" Loading epsilon from file [ "<<_EPSILON_FN<<" ] ..."<<endl;
	string tempstr;
	vector<string> tempvec;
	int EP_SIG_X=0;
	int EP_SIG_Y=0;
	int res_num=0;
	for(EP_SIG_X=0; EP_SIG_X<_NUM_chains; EP_SIG_X++){
		res_num+=_system_.chains[EP_SIG_X].nresidues;
	}
	while( getline(epsilon_ifstream, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		if(EP_SIG_Y==0) {
			EP_SIG_Y=tempvec.size();
		}
		if( tempvec.size()!=0 ) {
			EP_SIG_X+=1;
		}
	}
	epsilon_ifstream.close();
	if( EP_SIG_X<res_num ) {
		cout<<" _epsilon loading error: resnum="<<res_num<<", but N_vepsilone_x="<<EP_SIG_X<<";"<<endl;
		cout<<" check your [ "<<_EPSILON_FN<<" ];"<<endl;
		cout<<" exit~"<<endl;
		exit(LOGICERROR);
	} else {
		EP_SIG_X=res_num;
	}	
	if( EP_SIG_Y<res_num ) {
		cout<<" _epsilon loading error: resnum="<<res_num<<", but N_vepsilone_y="<<EP_SIG_Y<<";"<<endl;
		cout<<" check your [ "<<_EPSILON_FN<<" ];"<<endl;
		cout<<" exit~"<<endl;
		exit(LOGICERROR);
	} else {
		EP_SIG_Y=res_num;
	}

	CMyArray<double> EPSILON;
	CMyArray<double> SIGMA;
	EPSILON.Build(EP_SIG_X, EP_SIG_Y);
    SIGMA.Build(EP_SIG_X, EP_SIG_Y);

	ifstream epsilon_ifstream_1(_EPSILON_FN.c_str());
	int xx=0;
	int yy=0;
	while( getline(epsilon_ifstream_1, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		for(yy=0; yy<EP_SIG_Y; yy++ ) {
			EPSILON.pArray[xx][yy]=atof(tempvec[yy].c_str());
			cout<<setw(6)<<EPSILON.pArray[xx][yy]<<" ";
		}
		xx++;
		cout<<endl;
		if(xx>=EP_SIG_X) {
			break;
		}
	}
	epsilon_ifstream_1.close();
	
	ifstream sigma_ifstream(_SIGMA_FN.c_str());
	if(sigma_ifstream==NULL) {
		cout<<" Error: Can not open file: "<<sigma_ifstream<<endl;
		exit(IOERROR);
	}
	cout<<" Loading sigma from file [ "<<_SIGMA_FN<<" ] ..."<<endl;
	xx=0;
	yy=0;
	while( getline(sigma_ifstream, tempstr) ) {
		tempvec=Split(BFilter(tempstr));
		for(yy=0; yy<EP_SIG_Y; yy++ ) {
			SIGMA.pArray[xx][yy]=atof(tempvec[yy].c_str());
			cout<<setw(6)<<SIGMA.pArray[xx][yy]<<" ";
		}
		xx++;
		cout<<endl;
		if(xx>=EP_SIG_X) {
			break;
		}
	}
	sigma_ifstream.close();
	if(xx!=EP_SIG_X) {
		cout<<" when read _sigma.pls ... "<<endl;
		cout<<" xx = "<<xx<<endl;
		cout<<" EP_SIG_X = "<<EP_SIG_X<<endl;
		cout<<" error, check file and try again!"<<endl;
		exit(LOGICERROR);
	}
	/*yy=0;
	for(xx=0; xx<_NUM_chains; xx++){
		yy+=_system_.chains[xx].nresidues;
	}
	if( (EP_SIG_Y < yy) || (EP_SIG_X < yy) ) {
		cout<<" EP_SIG_X="<<EP_SIG_X<<endl;
		cout<<" EP_SIG_Y="<<EP_SIG_Y<<endl;
		cout<<" _system_.nresidues="<<yy<<endl;
		cout<<" error, check ur programs again!"<<endl;
	}*/
	CMyArray<double> Radius_cut;
	Radius_cut.Build(EP_SIG_X, EP_SIG_Y);
	cout<<" Loading r_cut from file [ _rcut.pls ] ... "<<endl;
	ifstream r_cut_istream("_rcut.pls");
	if(r_cut_istream==NULL) {
		cout<<" can not open [ _rcut.pls ], all set to infinity! "<<endl;
		for(xx=0; xx<EP_SIG_X; xx++ ) {
			for(yy=0; yy<EP_SIG_Y; yy++ ) {
				Radius_cut.pArray[xx][yy]=pow(_MAX_DOUBLE, 1.0/12.0);
				cout<<setw(6)<<Radius_cut.pArray[xx][yy]<<" ";
			}
			cout<<endl;
		}
	} else {
		xx=0;
		yy=0;
		while( getline(r_cut_istream, tempstr) ) {
			tempvec=Split(BFilter(tempstr));
			for(yy=0; yy<EP_SIG_Y; yy++ ) {
				if(tempvec[yy]==string("#")) {
					Radius_cut.pArray[xx][yy]=pow(2.0,1.0/6.0)*SIGMA.pArray[xx][yy];
				} else if(tempvec[yy]==string("@")) {
					Radius_cut.pArray[xx][yy]=pow(_MAX_DOUBLE, 1.0/12.0);
				} else {
					Radius_cut.pArray[xx][yy]=atof(tempvec[yy].c_str());
				}
				cout<<setw(6)<<Radius_cut.pArray[xx][yy]<<" ";
			}
			xx++;
			cout<<endl;
			if(xx>=EP_SIG_X) { break; }
		}
	}
	r_cut_istream.close();
	if(xx!=EP_SIG_X) {
		cout<<" when read _rcut.pls ... "<<endl;
		cout<<" xx = "<<xx<<endl;
		cout<<" EP_SIG_X = "<<EP_SIG_X<<endl;
		cout<<" error, check file and try again!"<<endl;
		exit(LOGICERROR);
	}
	int i=0;
	int j=0;
	int ep_sig_x=0;
	int ep_sig_y=0;
	for(i=1; i<_SIZE_memo; i++) {
		ep_sig_x=_INDEX_RES_ATM[i];
		for(j=1; j<_SIZE_memo; j++) {
			ep_sig_y=_INDEX_RES_ATM[j];
			_EPSILON_eachatom.pArray[i][j]=EPSILON.pArray[ep_sig_x][ep_sig_y];
			_SIGMA_eachatom.pArray[i][j]=SIGMA.pArray[ep_sig_x][ep_sig_y];
			_SIGMA3_eachatom.pArray[i][j]=TRIPLE(_SIGMA_eachatom.pArray[i][j]);
			_SIGMA6_eachatom.pArray[i][j]=_SIGMA3_eachatom.pArray[i][j]*_SIGMA3_eachatom.pArray[i][j];
			_SIGMA9_eachatom.pArray[i][j]=_SIGMA3_eachatom.pArray[i][j]*_SIGMA6_eachatom.pArray[i][j];
			_SIGMA12_eachatom.pArray[i][j]=_SIGMA6_eachatom.pArray[i][j]*_SIGMA6_eachatom.pArray[i][j];
			_R_cut_RR_eachatom.pArray[i][j]=Radius_cut.pArray[ep_sig_x][ep_sig_y]*Radius_cut.pArray[ep_sig_x][ep_sig_y];
			_DIS2=_R_cut_RR_eachatom.pArray[i][j];
			_DIS6=TRIPLE(_R_cut_RR_eachatom.pArray[i][j]);
			_DIS12=_DIS6*_DIS6;
			
			_E_cut_RR_eachatom.pArray[i][j]=Energy_LJ(_EPSILON_eachatom.pArray[i][j],_SIGMA12_eachatom.pArray[i][j],_SIGMA6_eachatom.pArray[i][j],_DIS12,_DIS6);
			
			    //cout<<i<<" "<<j<<" succ!"<<endl;
			/*if( _system_.atoms[i-1].residuename==string("ROD") && _system_.atoms[j-1].residuename==string("ROD") )
			{
				_IF_RR_eachatom.pArray[i][j]=1;//rod is a bond, not impossible a single atom;
			}
			else
			{
				_IF_RR_eachatom.pArray[i][j]=0;//rod is a bond, not impossible a single atom;
			}*/
			/*cout<<setw(4)<<_R_cut_RR_eachatom.pArray[i][j]<<" "<<_E_cut_RR_eachatom.pArray[i][j]
				<<" "<<_IF_RR_eachatom.pArray[i][j]<<endl;*/
		}
		//cout<<endl;
	}


	double sigma_max;
	int num_inteval=1000;
	double dis_inteval;
	double t_dis;
	double t_dis_2;
	double t_dis_6;
	double t_dis_12;
	double t_ener;
	double temp_sig6,temp_sig12;
	double temp_r6,temp_r12;
	int t_res=0;
	string i_name;
	string j_name;
	//int i=0;
	//int j=0;
	cout<<" [ ---------- gnuploting the potential profile ---------- ]"<<endl;
	ofstream pararecord("pararecord.par");
	for(ep_sig_x=0; ep_sig_x<EP_SIG_X; ep_sig_x++) {
		for(ep_sig_y=ep_sig_x; ep_sig_y<EP_SIG_Y; ep_sig_y++) {
			t_res=0;
			for(i=0;i<_NUM_chains;i++) {
				for(j=0;j<_system_.chains[i].nresidues;j++) {
					if( (t_res+j) == ep_sig_x ) {
						i_name=_system_.chains[i].chainname+string("_")+_system_.chains[i].residues[j].residuename;
						break;
					}
				}
				t_res+=_system_.nchains;
			}
			t_res=0;
			for(i=0;i<_NUM_chains;i++) {
				for(j=0;j<_system_.chains[i].nresidues;j++) {
					if( (t_res+j) == ep_sig_y ) {
						j_name=_system_.chains[i].chainname+string("_")+_system_.chains[i].residues[j].residuename;
						break;
					}
				}
				t_res+=_system_.nchains;
			}
			i_name=i_name+string("_")+j_name;
			ofstream para4gnuplot("para4gnuplot.par");
			para4gnuplot<<"fn=\'"<<i_name<<"\'\n"
			            <<"epsilon="<<showpoint<<setprecision(22)<<EPSILON.pArray[ep_sig_x][ep_sig_y]<<"\n"
			            <<"sigma="<<showpoint<<setprecision(22)<<SIGMA.pArray[ep_sig_x][ep_sig_y]<<"\n"
			            <<"cut="<<showpoint<<setprecision(22)<<Radius_cut.pArray[ep_sig_x][ep_sig_y]<<endl;
			pararecord<<"fn=\'"<<i_name<<"\'\n"
			          <<"epsilon="<<showpoint<<setprecision(22)<<EPSILON.pArray[ep_sig_x][ep_sig_y]<<"\n"
			          <<"sigma="<<showpoint<<setprecision(22)<<SIGMA.pArray[ep_sig_x][ep_sig_y]<<"\n"
			          <<"cut="<<showpoint<<setprecision(22)<<Radius_cut.pArray[ep_sig_x][ep_sig_y]<<"\n"<<"\n";
			para4gnuplot.close();
			i_name=i_name+string(".nrg");
			ofstream outener( i_name.c_str() );
			sigma_max=4.0*SIGMA.pArray[ep_sig_x][ep_sig_y];
			dis_inteval=sigma_max/num_inteval;
			temp_sig6=HEXIAL(SIGMA.pArray[ep_sig_x][ep_sig_y]);
			temp_sig12=temp_sig6*temp_sig6;
			temp_r6=HEXIAL(Radius_cut.pArray[ep_sig_x][ep_sig_y]);
			temp_r12=temp_r6*temp_r6;
			for(xx=1; xx<=num_inteval; xx++) {
				t_dis=xx*dis_inteval;
				t_dis_2=t_dis*t_dis;
				t_dis_6=TRIPLE(t_dis_2);
				t_dis_12=t_dis_6*t_dis_6;
				t_ener=0.0;
				if(t_dis<Radius_cut.pArray[ep_sig_x][ep_sig_y]) {
					/*if(ep_sig_x==1 && ep_sig_y==1 && t_dis>Radius_cut.pArray[ep_sig_x][ep_sig_y]-2*dis_inteval) {
						cout<<Radius_cut.pArray[ep_sig_x][ep_sig_y]<<" "<<pow(2,1/6.0)<<" "
						    <<EPSILON.pArray[ep_sig_x][ep_sig_y]<<" "
						    <<SIGMA.pArray[ep_sig_x][ep_sig_y]<<" "
						    <<temp_sig<<" "<<SQUARE(temp_sig)<<" "
						    <<temp_r<<" "<<SQUARE(temp_r)<<" "<<2*( pow( 1/pow(2,1/6.0), 12 ) - pow( 1/pow(2,1/6.0), 6) )<<" "
						    <<Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],SQUARE(temp_sig),temp_sig,t_dis_12,t_dis_6)<<" "
						    <<Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],SQUARE(temp_sig),temp_sig,SQUARE(temp_r),temp_r)<<" "
						    <<Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],SQUARE(temp_sig),temp_sig,t_dis_12,t_dis_6)-Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],SQUARE(temp_sig),temp_sig,SQUARE(temp_r),temp_r)<<endl;
						break;
					}*/
					t_ener=Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],temp_sig12,temp_sig6,t_dis_12,t_dis_6)
				    	  -Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],temp_sig12,temp_sig6,temp_r12,temp_r6);
					/*if(ep_sig_x != ep_sig_y) {
						cout<<t_dis<<" "<<Radius_cut.pArray[ep_sig_x][ep_sig_x]<<" "<<Energy_LJ(EPSILON.pArray[ep_sig_x][ep_sig_y],temp_sig12,temp_sig6,temp_r12,temp_r6)<<" "<<t_ener<<endl;
					}*/
				}
				outener<<t_dis<<"\t"<<t_ener<<"\n";
			}
			outener.close();
			system("gnuplot < draw.gpl");
		}
	}
	pararecord.close();
	cout<<" [ ---------- gnuploting the potential profile ---------- ]"<<endl;
	cout<<" [ NOTE ]: you can check your energy profile now... those .eps files. "<<endl;
	cout<<" "<<endl;
	
	/*vector<int> _SEQ_agl_rod_atoms_temp;
	vector<int> _SEQ_agl_rodpairs_temp;
	for(i=1; i<tempSize; i++)
	{		
		if( _TYPE_atom_all_atom[i]==0 && _IF_RR_eachatom.pArray[i][i-1]==1 && _IF_RR_eachatom.pArray[i][i+1]==1 )
		{
			_IF_RR_calc_or_not[i]=1;
			_Cos_RR_eachatom[i]=cos_angle( _XX[i-1], _YY[i-1], _ZZ[i-1], 
			                               _XX[i],   _YY[i],   _ZZ[i],
								           _XX[i+1], _YY[i+1], _ZZ[i+1] );
			_SEQ_agl_rod_atoms_temp.push_back(i);
		}
		else
		{
			_IF_RR_calc_or_not[i]=0;
			if( _TYPE_atom_all_atom[i]==0 )
			{
				if( _IF_RR_eachatom.pArray[i][i-1]==0 && _IF_RR_eachatom.pArray[i][i+1]==1 )
				{
					_SEQ_agl_rodpairs_temp.push_back(i);
				}
				else if(_IF_RR_eachatom.pArray[i][i-1]==1 && _IF_RR_eachatom.pArray[i][i+1]==0)
				{
					_SEQ_agl_rodpairs_temp.push_back(i);
				}
				else
				{
				}
				//cout<<_IF_RR_eachatom.pArray[i][i-1]<<" "<<_IF_RR_eachatom.pArray[i][i+1]<<" ";
			}
			else if( _TYPE_atom_all_atom[i]==-1 )
			{
				if( _IF_RR_eachatom.pArray[i][i+1]==1 )//rod is a bond, not impossible a single atom;
				{
					_SEQ_agl_rodpairs_temp.push_back(i);
				}
				else
				{
				}
				//cout<<_IF_RR_eachatom.pArray[i][i+1]<<" ";
			}
			else if(_TYPE_atom_all_atom[i]==1)
			{
				if( _IF_RR_eachatom.pArray[i][i-1]==1 )//rod is a bond, not impossible a single atom;
				{
					_SEQ_agl_rodpairs_temp.push_back(i);
				}
				else
				{
				}
				//cout<<_IF_RR_eachatom.pArray[i][i-1]<<" ";
			}
			//cout<<i<<" ch: "<<_SEQ_agl_rodpairs_temp[_SEQ_agl_rodpairs_temp.size()-1]<<endl;
		}
	}
	_NUM_agl_rod_atoms=_SEQ_agl_rod_atoms_temp.size();
	if(_NUM_agl_rod_atoms!=0)
	{
		_SEQ_agl_rod_atoms=new int[_NUM_agl_rod_atoms];
		MEMOSETZERO(_SEQ_agl_rod_atoms, sizeof(int)*_NUM_agl_rod_atoms);
		for(i=0; i<_NUM_agl_rod_atoms; i++)
		{
			_SEQ_agl_rod_atoms[i]=_SEQ_agl_rod_atoms_temp[i];
			//cout<<setw(3)<<" "<<_SEQ_agl_rod_atoms[i];
		}
		//cout<<endl;
		_SEQ_agl_rod_atoms_temp.clear();
	}
	//////////////////////////////////////////
	if(_SEQ_agl_rodpairs_temp.size()%2!=0)
	{
		cout<<" _SEQ_agl_rodpairs_temp.size() = "<<_SEQ_agl_rodpairs_temp.size()<<endl;
		cout<<" error, check program and try again!"<<endl;
		exit(LOGICERROR);
	}
	_NUM_agl_rodresidues=_SEQ_agl_rodpairs_temp.size()/2;
	_NUM_agl_rodpairs=_NUM_agl_rodresidues*(_NUM_agl_rodresidues-1)/2;
	if(_NUM_agl_rodpairs!=0)
	{
		_Cos_RR_eachrodpair=new double[_NUM_agl_rodpairs];
	}
	if(_NUM_agl_rodresidues!=0)
	{
		_SEQ_agl_rodpairs=new int[_NUM_agl_rodresidues*2];
		for(i=0; i<_NUM_agl_rodresidues; i++)
		{
			_SEQ_agl_rodpairs[2*i]=_SEQ_agl_rodpairs_temp[2*i];
			_SEQ_agl_rodpairs[2*i+1]=_SEQ_agl_rodpairs_temp[2*i+1];
			cout<<" rod["<<setw(2)<<i<<"]"
				<<": s[ "<<setw(3)<<_SEQ_agl_rodpairs[2*i]<<" <---> "
				<<setw(3)<<_SEQ_agl_rodpairs[2*i+1]<<" ]e "<<endl;
		}
		if(_NUM_agl_rodpairs!=0)
		{
			int num=0;
			for(i=0; i<_NUM_agl_rodresidues-1; i++)
			{
				for(j=i+1; j<_NUM_agl_rodresidues; j++)
				{
					_Cos_RR_eachrodpair[num++]=cos_angle( 
						_XX[_SEQ_agl_rodpairs[2*i]], _YY[_SEQ_agl_rodpairs[2*i]], _ZZ[_SEQ_agl_rodpairs[2*i]], 
						_XX[_SEQ_agl_rodpairs[2*i+1]],   _YY[_SEQ_agl_rodpairs[2*i+1]],   _ZZ[_SEQ_agl_rodpairs[2*i+1]],
						_XX[_SEQ_agl_rodpairs[2*j]], _YY[_SEQ_agl_rodpairs[2*j]], _ZZ[_SEQ_agl_rodpairs[2*j]],
						_XX[_SEQ_agl_rodpairs[2*j+1]], _YY[_SEQ_agl_rodpairs[2*j+1]], _ZZ[_SEQ_agl_rodpairs[2*j+1]]);
				}
			}
			if(num!=_NUM_agl_rodpairs)
			{
				cout<<" num!=_NUM_agl_rodpairs "<<endl;
				cout<<" error, check program and try again!"<<endl;
				exit(LOGICERROR);
			}
		}
		_SEQ_agl_rodpairs_temp.clear();
	}*/
	//double
	EPSILON.Release();
	SIGMA.Release();
	Radius_cut.Release();
}
///////////////////////////
///////////////////////////
bool cmc::make_change() {
	if( _INDEX_chn_or_not ) {// if atom belongs to a chain;
		Omega=_PI_double*rand_seed(iseed_angle);
		sin_Omega=sin(Omega);
		cos_Omega=cos(Omega);
		if(!_BF_flag) { // without bond fluctuation
			if( _TYPE_atom_ind == 0 ) {// if atom is in the middle of the chain;  	
				//vector AB -> ;  
				dc[0]=_XX[_INDEX_chosen+1]-_XX[_INDEX_chosen-1];
				dc[1]=_YY[_INDEX_chosen+1]-_YY[_INDEX_chosen-1];
				dc[2]=_ZZ[_INDEX_chosen+1]-_ZZ[_INDEX_chosen-1];
				dd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];

				ec[0]=_XX[_INDEX_chosen]-_XX[_INDEX_chosen-1];
				ec[1]=_YY[_INDEX_chosen]-_YY[_INDEX_chosen-1];
				ec[2]=_ZZ[_INDEX_chosen]-_ZZ[_INDEX_chosen-1];
				ee=ec[0]*ec[0]+ec[1]*ec[1]+ec[2]*ec[2];

				fc[0]=_XX[_INDEX_chosen+1]-_XX[_INDEX_chosen];
				fc[1]=_YY[_INDEX_chosen+1]-_YY[_INDEX_chosen];
				fc[2]=_ZZ[_INDEX_chosen+1]-_ZZ[_INDEX_chosen];
				ff=fc[0]*fc[0]+fc[1]*fc[1]+fc[2]*fc[2];

				//just for safe, although this will be impossible in a real system.
				//printf(" dd=%32.30f\n", dd);
				if(dd<1e-12) {
					/*cout<<"TempEner="<<TempEner<<endl;
					cout<<"_ENER_total="<<_ENER_total<<endl;
					cout<<"tempindex_judge="<<tempindex_judge<<endl;
					cout<<"tempindex_judge_bak="<<tempindex_judge_bak<<endl;*/
					factor=0.5;
					cout<<" factor="<<factor<<" truly happened!"<<endl;
				} else {
					factor=(ee-ff)/2.0/dd+0.5;
				}

				//vector AD -> ;
				uc=dc[0]*factor;//+_XX[_INDEX_chosen-1];
				vc=dc[1]*factor;//+_YY[_INDEX_chosen-1];
				wc=dc[2]*factor;//+_ZZ[_INDEX_chosen-1];
				
				//vector DC -> ; DC = AC - AD ;  
				vecDC[0]=ec[0]-uc;
				vecDC[1]=ec[1]-vc;
				vecDC[2]=ec[2]-wc;

				dd_sqrt=sin_Omega/sqrt(dd); //not sqrt(d) but for the convenience for next step calc;

				//coorRES=Add_p(C_p(cos_Omega,vecDC), C_p(sin_Omega/sqrt(dd), p_X_p(dc, vecDC)));
				_XX[_INDEX_chosen] = _XX[_INDEX_chosen-1] + uc + vecDC[0]*cos_Omega + (dc[1]*vecDC[2]-dc[2]*vecDC[1]) * dd_sqrt;
				_YY[_INDEX_chosen] = _YY[_INDEX_chosen-1] + vc + vecDC[1]*cos_Omega + (dc[2]*vecDC[0]-dc[0]*vecDC[2]) * dd_sqrt;
				_ZZ[_INDEX_chosen] = _ZZ[_INDEX_chosen-1] + wc + vecDC[2]*cos_Omega + (dc[0]*vecDC[1]-dc[1]*vecDC[0]) * dd_sqrt;
				
		
			} else {// if atom is at the front or at the end of the chain; 
				THETA=asin(rand_seed(iseed_angle)*2.0-1.0)+_PI_half;
				//cout<<" THE="<<THE<<endl;
				sin_THETA=sin(THETA);
				cos_THETA=cos(THETA);

				tc[0]=_XX[_INDEX_chosen-_TYPE_atom_ind]-_XX[_INDEX_chosen];
				tc[1]=_YY[_INDEX_chosen-_TYPE_atom_ind]-_YY[_INDEX_chosen];
				tc[2]=_ZZ[_INDEX_chosen-_TYPE_atom_ind]-_ZZ[_INDEX_chosen];

				tl=sqrt(tc[0]*tc[0]+tc[1]*tc[1]+tc[2]*tc[2]);
				//cout<<" tl="<<tl<<endl;
				sin_THETA=sin_THETA*tl;
				_XX[_INDEX_chosen]=_XX[_INDEX_chosen-_TYPE_atom_ind]+cos_Omega*sin_THETA;
				_YY[_INDEX_chosen]=_YY[_INDEX_chosen-_TYPE_atom_ind]+sin_Omega*sin_THETA;
				_ZZ[_INDEX_chosen]=_ZZ[_INDEX_chosen-_TYPE_atom_ind]+cos_THETA*tl;
			}
		} else {//with bond fluctuation
			if( _TYPE_atom_ind == 0 ) {// if atom is in the middle of the chain;  			
				dc[0]=_XX[_INDEX_chosen+1]-_XX[_INDEX_chosen-1];
				dc[1]=_YY[_INDEX_chosen+1]-_YY[_INDEX_chosen-1];
				dc[2]=_ZZ[_INDEX_chosen+1]-_ZZ[_INDEX_chosen-1];
				dd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];
				dd_sqrt=sqrt(dd);
								
				ee=rmin+rand_seed(iseed_len1)*ddelta;
				r2max=min((dd_sqrt+ee), rmax);
				r2min=max(fabs(dd_sqrt-ee), rmin);
				if(r2max>r2min) {
					//cout<<"r2max="<<r2max<<"::r2min="<<r2min<<endl;
					ff=r2min+rand_seed(iseed_len1)*(r2max-r2min);
					//cout<<ff<<endl;
				} else { return false; }				
				
				/*if( ee<rmin || ff<rmin )//for check!
				{
					printf("rmax=%e, rmin=%e\n", rmax, rmin);
					printf("ee=%e\n", ee);
					printf("sqrt(dd)=%e\n", dd_sqrt);
					printf("r2max=%e, r2min=%e\n", r2max, r2min);
					printf("ff=%e\n", ff);	
				}*/
					
				ee=ee*ee;
				ff=ff*ff;
				
				//just for safe, although this will be impossible in a real system.
				//printf(" dd=%32.30f\n", dd);
				if(dd<1e-12)
				{
					factor=0.5;
					cout<<" factor="<<factor<<" truly happened!"<<endl;
				}
				else
				{
					factor=(ee-ff+dd)/2.0/dd;
				}
				tl=sqrt(ee-dd*factor*factor);
				if(tl!=tl)//when e+f=d; tl=nan; this is important
				{
					tl=0.0;
				}

				uc=dc[0]*factor+_XX[_INDEX_chosen-1];
				vc=dc[1]*factor+_YY[_INDEX_chosen-1];
				wc=dc[2]*factor+_ZZ[_INDEX_chosen-1];
					
				TempBX=_XX[_INDEX_chosen+1]-uc;
				TempBY=_YY[_INDEX_chosen+1]-vc;
				TempBZ=_ZZ[_INDEX_chosen+1]-wc;

				Len_CB_xy_2=TempBX*TempBX+TempBY*TempBY;
				Len_CB_xy=sqrt(Len_CB_xy_2);
				Len_CB=sqrt(Len_CB_xy_2+TempBZ*TempBZ);

				if( fabs(Len_CB_xy)<1e-12 )
				{
					sin_PHI=1.0;
					cos_PHI=0.0;
					if( fabs(Len_CB)<1e-12 )
					{
						sin_THETA=1.0;
						cos_THETA=0.0;
					}
					else
					{						
						sin_THETA=Len_CB_xy/Len_CB;
						cos_THETA=TempBZ/Len_CB;
					}
				}
				else
				{
					sin_PHI=TempBY/Len_CB_xy;
					cos_PHI=TempBX/Len_CB_xy;
					sin_THETA=Len_CB_xy/Len_CB;
					cos_THETA=TempBZ/Len_CB;
				}
				
				_XX[_INDEX_chosen]=uc+( cos_PHI*cos_THETA*cos_Omega - sin_PHI*sin_Omega )*tl;
				_YY[_INDEX_chosen]=vc+( sin_PHI*cos_THETA*cos_Omega + cos_PHI*sin_Omega )*tl;
				_ZZ[_INDEX_chosen]=wc-( sin_THETA*cos_Omega )*tl;
			}
			if( _TYPE_atom_ind == -1 || _TYPE_atom_ind == 1 ) {// if atom is at the front or the end of the chain; 
				cos_THETA=1.0-rand_seed(iseed_angle)*2.0; //wrong
				sin_THETA=sqrt(1-cos_THETA*cos_THETA);

				tl=rmin+rand_seed(iseed_len1)*ddelta;

				_XX[_INDEX_chosen]=_XX[_INDEX_chosen-_TYPE_atom_ind]+cos_Omega*sin_THETA*tl;
				_YY[_INDEX_chosen]=_YY[_INDEX_chosen-_TYPE_atom_ind]+sin_Omega*sin_THETA*tl;
				_ZZ[_INDEX_chosen]=_ZZ[_INDEX_chosen-_TYPE_atom_ind]+cos_THETA*tl;
			}
		}
	} else {// if atom does not belongs to a chain;	
	}
	_MC_NUM_TOT+=1; /*************  for check program ************/
	return true;
}
///////////////////////////////////////////////////////////////
///////////////////////////
bool cmc::chck_bond_len()//after each make mapping!
{
	if(!_Enery_initialized) {
		cout<<" energy not initialized, can not check bond length..."<<endl;
		return LOGICERROR;
	}
	int i=0;
	double neighbor_len=0.0;
	for(i=1; i<_NUM_atoms; i++) {
		if(_INDEX_chn_atom[i+1]==_INDEX_chn_atom[i]) {
			neighbor_len=sqrt(_DIS2_eachatom.pArray[i][i+1]);
			if( !_BF_flag ) {
				if( fabs(neighbor_len-_BOND_length)>1e-3 ) {
					cout<<endl<<" distance between "<<i<<" and "<<i+1<<" atoms has some problem: "<<endl;
					printf(" neighbor_len=%2.18lf\n", neighbor_len);
					printf(" _BOND_length=%2.18lf\n", _BOND_length);
					printf(" neighbor_len-_BOND_length=%e\n", neighbor_len-_BOND_length);
					cout<<" error: |neighbor_len-_BOND_length|>1e-6, check your program."<<endl;
					return false;
				}
			}
			else if( fabs(neighbor_len-_BOND_length)>_BOND_delta  ) {
				cout<<endl<<" distance between "<<i<<" and "<<i+1<<" atoms has some problem: "<<endl;
				printf(" neighbor_len=%e\n", neighbor_len);
				printf(" _BOND_length=%e\n", _BOND_length);
				printf(" neighbor_len-_BOND_length=%e\n", neighbor_len-_BOND_length);
				cout<<" error: |neighbor_len-_BOND_length|>_BOND_delta, check your program."<<endl;
				return false;
			}
		}
	}
	return true;
}
///////////////////////////
void cmc::make_choice(const int INDEX_coor) {	
	if( (INDEX_coor>_NUM_atoms) || (INDEX_coor<=0) ) {
		cout<<"Error: INDEX_coor="<<INDEX_coor<<", should be in [1, "<<_NUM_atoms<<"]"<<endl;
		exit(LOGICERROR);
	}
	_INDEX_chosen           = INDEX_coor;
	_INDEX_chn_ind          =_INDEX_chn_atom[_INDEX_chosen];        // for the chosen one;
	_SIZE_of_chosen_chn     =_SIZE_of_chn[_INDEX_chn_ind];          // for the chosen one; * different!
	//_INDEX_chnhead          =_INDEX_CHN_HEAD[_INDEX_chn_ind];
	//_INDEX_chntail          =_INDEX_CHN_TAIL[_INDEX_chn_ind];
	_INDEX_chn_or_not       =_INDEX_chn_or_not_atom[_INDEX_chosen]; // for the chosen one;
	//_INDEX_atm_in_chn_ind   =_INDEX_atm_in_chn[_INDEX_chosen]; // for the chosen one;
	_TYPE_atom_ind          =_TYPE_atom[_INDEX_chosen];        // for the chosen one;
	_Index_lneighbor_ind    =_INDEX_LNEIGHBOR[_INDEX_chosen];
	_Index_rneighbor_ind    =_INDEX_RNEIGHBOR[_INDEX_chosen]+1;
	//// [preparation for change...]
	_XX[0]=_XX[_INDEX_chosen];
	_YY[0]=_YY[_INDEX_chosen];
	_ZZ[0]=_ZZ[_INDEX_chosen];
	//calc_energy_chosen_atom();// not important any more in this detailed energy edition !!!!!!!!!!!!!!!!!
	for(int i=1; i<_SIZE_memo; i++) {
		_ENER_LJ_chosenatom_backup[i]=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i];
		_DIS2_chosenatom_backup[i]=_DIS2_eachatom.pArray[_INDEX_chosen][i];
		_DIS_x_chosenatom_backup[i]=_DIS_x_eachatom.pArray[_INDEX_chosen][i];
		_DIS_y_chosenatom_backup[i]=_DIS_y_eachatom.pArray[_INDEX_chosen][i];
		_DIS_z_chosenatom_backup[i]=_DIS_z_eachatom.pArray[_INDEX_chosen][i];
	}	
	_ENER_BF_chosenatom_backup=_ENER_BF_eachatom[_INDEX_chosen];
}
///////////////////////////
void cmc::make_judge() {	
    calc_energy_chosen_atom();
	_ENER_delta=0.0;
	for(int i=1; i<_SIZE_memo; i++) {
		_ENER_delta+=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i]-_ENER_LJ_chosenatom_backup[i]
		            +_ENER_BF_eachatom[_INDEX_chosen]-_ENER_BF_chosenatom_backup;
	}
	/*if( (_INDEX_chosen==1) || (_INDEX_chosen==_NUM_atoms) )
	{
		_ENER_delta+=_ENER_ele_new-_ENER_ele_old;
	}*/
	if(_ENER_delta<=0.0) {
		tempnumer_ret=2;//succ!
	} else {
		if( exp(-_ENER_delta*_TEMPERATURE_REP) > rand_seed(iseed_rand) ) {
			tempnumer_ret=2;// succ;
		} else {
			tempnumer_ret=1;// fail;
		}
	}
	//}
}
///////////////////////////
void cmc::make_accept() {
	_ENER_total+=_ENER_delta;
	/*if( (_INDEX_chosen==1) || (_INDEX_chosen==_NUM_atoms) )
	{
		_ENER_ele_old=_ENER_ele_new;
		//cout<<_ENER_ele_old/_ENER_total<<endl;//for test!
	}*/
	for(int i=1; i<_SIZE_memo; i++) {
		_ENER_LJ_eachatom.pArray[i][_INDEX_chosen]=_ENER_LJ_eachatom.pArray[_INDEX_chosen][i];
		_DIS2_eachatom.pArray[i][_INDEX_chosen]=_DIS2_eachatom.pArray[_INDEX_chosen][i];
		_DIS_x_eachatom.pArray[i][_INDEX_chosen]=_DIS_x_eachatom.pArray[_INDEX_chosen][i];
		_DIS_y_eachatom.pArray[i][_INDEX_chosen]=_DIS_y_eachatom.pArray[_INDEX_chosen][i];
		_DIS_z_eachatom.pArray[i][_INDEX_chosen]=_DIS_z_eachatom.pArray[_INDEX_chosen][i];
	}
	if(_FLAG_rg2) {
		_COM_x[_INDEX_chn_ind]+=_XX[_INDEX_chosen]-_XX[0];
		_COM_y[_INDEX_chn_ind]+=_YY[_INDEX_chosen]-_YY[0];
		_COM_z[_INDEX_chn_ind]+=_ZZ[_INDEX_chosen]-_ZZ[0];
	}
	//make_mapping_individual();//if just _XX when calc_ener, useless!
	_MC_NUM_SUC+=1;/*************  for check program ************/
}
///////////////////////////
void cmc::make_reject_all()
{
	_XX[_INDEX_chosen]=_XX[0];
	_YY[_INDEX_chosen]=_YY[0];
	_ZZ[_INDEX_chosen]=_ZZ[0];
	/*if( (_INDEX_chosen==1) || (_INDEX_chosen==_NUM_atoms) )
	{
		_ENER_ele_new=_ENER_ele_old;
	}*/
	for(int i=1; i<_SIZE_memo; i++) {
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][i]=_ENER_LJ_chosenatom_backup[i];
		_DIS2_eachatom.pArray[_INDEX_chosen][i]=_DIS2_chosenatom_backup[i];
		_DIS_x_eachatom.pArray[_INDEX_chosen][i]=_DIS_x_chosenatom_backup[i];
		_DIS_y_eachatom.pArray[_INDEX_chosen][i]=_DIS_y_chosenatom_backup[i];
		_DIS_z_eachatom.pArray[_INDEX_chosen][i]=_DIS_z_chosenatom_backup[i];
	}
	_ENER_BF_eachatom[_INDEX_chosen]=_ENER_BF_chosenatom_backup;
	_MC_NUM_FIL+=1;/*************  for check program ************/
}
///////////////////////////
void cmc::make_reject_only_move() {
	_XX[_INDEX_chosen]=_XX[0];
	_YY[_INDEX_chosen]=_YY[0];
	_ZZ[_INDEX_chosen]=_ZZ[0];
	//_MC_NUM_FIL+=1;/*************  for check program ************///no need!
}
///////////////////////////
int cmc::make_mcmove(const int INDEX_coor) {
	make_choice(INDEX_coor);
	if( make_change() ) {
		make_judge();//in judge _ENER_total can not be changed, should be in accept! 
		if(tempnumer_ret==2) {//succ==2
			//cout<<INDEX_coor<<endl;
			make_accept();
		} else if(tempnumer_ret==1) {//fail==1
			make_reject_all();
		} else if(tempnumer_ret==-1) {//false==-1, cos' no new energy calc!
			make_reject_only_move();
		}
		return tempnumer_ret;
	} else {
		//do make sure that there is no change in coordinates,
		//or you must add make_reject_only_move here!
		return -1; //not count, cuz this is the calculation problem, 
		           //actually there is a solution but we did not find it; 
	}
	
}
///////////////////////////
void cmc::make_mapping() {
	int i=0;
	int j=0;
	int k=0;
	int numerator=0;
	int sz_res=0;
	int sz_atm=0;
	for(i=0; i<_NUM_chains; i++) {
		sz_res=_system_.chains[i].nresidues;
		for(j=0; j<sz_res; j++) {
			sz_atm=_system_.chains[i].residues[j].natoms;
			for(k=0; k<sz_atm; k++) {
				numerator++;			
				_system_.chains[i].residues[j].atoms[k].x_coordinate=_XX[numerator];
				_system_.chains[i].residues[j].atoms[k].y_coordinate=_YY[numerator];
				_system_.chains[i].residues[j].atoms[k].z_coordinate=_ZZ[numerator];
			}
		}
	}
	if(numerator!=_NUM_atoms) {
		cout<<" Total_Num="<<numerator<<endl;
		cout<<" _NUM_atoms="<<_NUM_atoms<<endl;
		ErrorMSG("Total_Num!=_NUM_atoms :: cmc::make_mapping");
		exit(SIZEERROR);
	}
}
///////////////////////////
void cmc::fout_conformation(const string FILENAME_conf, bool CHCK_BOND_LEN_OR_NOT) {
	translate2box();
	make_mapping();
	if(CHCK_BOND_LEN_OR_NOT) chck_bond_len();//after each make mapping!
	_system_.writepdbinfo(FILENAME_conf.c_str(), _IFVERBOSE);
}
/*//////////////////////////
inline void cmc::trajectory_rec() {
	
	_OUTfile=fopen( _FNoutput,"ab");    
    //fwrite(&_ENER_total,sizeof(double),1,_OUTfile);
    fwrite(&_XX_rec,sizeof(double),_NUM_rec,_OUTfile);
    fwrite(&_YY_rec,sizeof(double),_NUM_rec,_OUTfile);
    fwrite(&_ZZ_rec,sizeof(double),_NUM_rec,_OUTfile);
    fclose(_OUTfile);
}
////////////////////////////
inline void cmc::recording() {
	////// add your code here /////
	for(int i=0;i<_SIZE_memo;i++) {
		_XX_rec[_I_eachstep+i]=_XX[i];
		_YY_rec[_I_eachstep+i]=_YY[i];
		_ZZ_rec[_I_eachstep+i]=_ZZ[i];
	}
}*/
////////////////////////////
void cmc::calc_energy_chosen_atom() {
	int k=0;
	//int l=0;
	int x=0;
	int y=0;
	int z=0;
	//int it_nb_start=_INDEX_atm_in_chn_ind>_NEIGHBOR_lj?(_INDEX_chosen-_NEIGHBOR_lj):(_INDEX_chosen-_INDEX_atm_in_chn_ind);
	//int it_nb_end=(_INDEX_atm_in_chn_ind+1+_NEIGHBOR_lj)<_SIZE_of_chosen_chn?(_INDEX_chosen+_NEIGHBOR_lj+1):(_INDEX_chosen+_SIZE_of_chosen_chn-_INDEX_atm_in_chn_ind);

	//for test	
	/*int it_nb_start_bak=0;
	int it_nb_end_bak=0;
	if( _INDEX_atm_in_chn_ind >= _NEIGHBOR_lj ) {
		it_nb_start_bak=_INDEX_chosen-_NEIGHBOR_lj;
	} else {
		it_nb_start_bak=_INDEX_chosen-_INDEX_atm_in_chn_ind;
	}
	if( _INDEX_atm_in_chn_ind+1+_NEIGHBOR_lj <= _SIZE_of_chosen_chn ) {
		it_nb_end_bak=_INDEX_chosen+_NEIGHBOR_lj+1;
	} else {
		it_nb_end_bak=_INDEX_chosen+_SIZE_of_chosen_chn-_INDEX_atm_in_chn_ind;
	}
	if( it_nb_start!=it_nb_start_bak || it_nb_end!=it_nb_end_bak ) {
		cout<<it_nb_start<<" : "<<it_nb_start_bak<<endl;
		cout<<it_nb_end<<" : "<<it_nb_end_bak<<endl;
	}*/
	//test ends here;

	if( _INDEX_chn_or_not ) {
		if( ( _SIZE_of_chosen_chn > 1 ) && _BF_flag ) {//if chain && with Bond Fluctuation
			if( _TYPE_atom_ind == -1 ) { //the start(-1) or the end(1) of the chain
				tempind=_INDEX_chosen+1;
				_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]=_XX[_INDEX_chosen]-_XX[tempind];
				_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]=_YY[_INDEX_chosen]-_YY[tempind];
				_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]=_ZZ[_INDEX_chosen]-_ZZ[tempind];
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_DIS2_eachatom.pArray[_INDEX_chosen][tempind]=_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[_INDEX_chosen][tempind];
				_ENER_BF_eachatom[_INDEX_chosen]=Energy_BF(_PARA_K,_DIS2_eachatom.pArray[_INDEX_chosen][tempind]);
				//enumer++;
			} else if (_TYPE_atom_ind == 1) { 
				tempind=_INDEX_chosen-1;
				_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]=_XX[_INDEX_chosen]-_XX[tempind];
				_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]=_YY[_INDEX_chosen]-_YY[tempind];
				_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]=_ZZ[_INDEX_chosen]-_ZZ[tempind];
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_DIS2_eachatom.pArray[_INDEX_chosen][tempind]=_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[_INDEX_chosen][tempind];
				_ENER_BF_eachatom[tempind]=Energy_BF(_PARA_K,_DIS2_eachatom.pArray[_INDEX_chosen][tempind]);
		    } else {
				tempind=_INDEX_chosen+1;
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]=_XX[_INDEX_chosen]-_XX[tempind];
				_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]=_YY[_INDEX_chosen]-_YY[tempind];
				_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]=_ZZ[_INDEX_chosen]-_ZZ[tempind];
				_DIS2_eachatom.pArray[_INDEX_chosen][tempind]=_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[_INDEX_chosen][tempind];
				_ENER_BF_eachatom[_INDEX_chosen]=Energy_BF(_PARA_K,_DIS2_eachatom.pArray[_INDEX_chosen][tempind]);
				
				tempind=_INDEX_chosen-1;
				//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
				_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]=_XX[_INDEX_chosen]-_XX[tempind];
				_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]=_YY[_INDEX_chosen]-_YY[tempind];
				_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]=_ZZ[_INDEX_chosen]-_ZZ[tempind];
				_DIS2_eachatom.pArray[_INDEX_chosen][tempind]=_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_x_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_y_eachatom.pArray[_INDEX_chosen][tempind]
				                                             +_DIS_z_eachatom.pArray[_INDEX_chosen][tempind]*_DIS_z_eachatom.pArray[_INDEX_chosen][tempind];
				_ENER_BF_eachatom[tempind]=Energy_BF(_PARA_K,_DIS2_eachatom.pArray[_INDEX_chosen][tempind]);
				
				//enumer++;
			}
		}
		//cout<<_INDEX_chosen<<"::"<<tempind<<endl;
		if ( ( _SIZE_of_chosen_chn > 2 ) && _AG_flag ) {//if chain and with angle calculation;
		}
		if ( ( _SIZE_of_chosen_chn > 3 ) && _DH_flag ) {//if chain and with dihedral calculation;	
		}
	} else { // if not chain, L-J potential calculation between "neighboring atoms";
		if( _SIZE_of_chosen_chn > 1 ) { // 2 atoms or more;
			for(k=_Index_lneighbor_ind; k<_INDEX_chosen; k++) {
				_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[_INDEX_chosen]-_XX[k]);
				_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[_INDEX_chosen]-_YY[k]);
				_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[_INDEX_chosen]-_ZZ[k]);
				_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                                +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                                +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
				for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
					for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
						for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
							_DIS_x=_PBL_X*x;
							_DIS_y=_PBL_Y*y;
							_DIS_z=_PBL_Z*z;
							_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
							                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
							                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
							                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
							_DIS6=TRIPLE(_DIS2);
							_DIS12=_DIS6*_DIS6; 
							if(_DIS2<_R_cut_RR_eachatom.pArray[_INDEX_chosen][k]) {
								_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachatom.pArray[_INDEX_chosen][k],_SIGMA12_eachatom.pArray[_INDEX_chosen][k],_SIGMA6_eachatom.pArray[_INDEX_chosen][k],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[_INDEX_chosen][k];
							}
						}
					}
				}
			}
			for(k=_INDEX_chosen+1; k<_Index_rneighbor_ind; k++) {
				_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[_INDEX_chosen]-_XX[k]);
				_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[_INDEX_chosen]-_YY[k]);
				_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[_INDEX_chosen]-_ZZ[k]);
				_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
				                                +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
				                                +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];
				_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
				for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
					for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
						for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
							_DIS_x=_PBL_X*x;
							_DIS_y=_PBL_Y*y;
							_DIS_z=_PBL_Z*z;
							_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
							                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
							                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
							                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
							_DIS6=TRIPLE(_DIS2);
							_DIS12=_DIS6*_DIS6;
							if(_DIS2<_R_cut_RR_eachatom.pArray[_INDEX_chosen][k]) {
								_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachatom.pArray[_INDEX_chosen][k],_SIGMA12_eachatom.pArray[_INDEX_chosen][k],_SIGMA6_eachatom.pArray[_INDEX_chosen][k],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[_INDEX_chosen][k];
							}
						}
					}
				}
			}
		}
	}

    //Lennard-Jones Potential Part; for intra chain ( not neighboring atoms ) and inter chain;
	//cout<<" [k:";
	for(k=1; k<_Index_lneighbor_ind; k++) {// must be '<', attention! 
		//cout<<" "<<k;
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[_INDEX_chosen]-_XX[k]);
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[_INDEX_chosen]-_YY[k]);
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[_INDEX_chosen]-_ZZ[k]);
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
		                                +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
		                                +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
		for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
			for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
				for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
					_DIS_x=_PBL_X*x;
					_DIS_y=_PBL_Y*y;
					_DIS_z=_PBL_Z*z;
					_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
					                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
					_DIS6=TRIPLE(_DIS2);
					_DIS12=_DIS6*_DIS6;
					if(_DIS2<_R_cut_RR_eachatom.pArray[_INDEX_chosen][k]) {
						_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachatom.pArray[_INDEX_chosen][k],_SIGMA12_eachatom.pArray[_INDEX_chosen][k],_SIGMA6_eachatom.pArray[_INDEX_chosen][k],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[_INDEX_chosen][k];
					}
				}
			}
		}
	}
	//cout<<" ("<<_INDEX_chosen<<")";
	for(k=_Index_rneighbor_ind; k<_SIZE_memo; k++) {// must be '<', attention! 
		//cout<<" "<<k;
		_DIS_x_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_X(_XX[_INDEX_chosen]-_XX[k]);
		_DIS_y_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Y(_YY[_INDEX_chosen]-_YY[k]);
		_DIS_z_eachatom.pArray[_INDEX_chosen][k]=DIS_PBC_Z(_ZZ[_INDEX_chosen]-_ZZ[k]);
		_DIS2_eachatom.pArray[_INDEX_chosen][k]=_DIS_x_eachatom.pArray[_INDEX_chosen][k]*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
		                                +_DIS_y_eachatom.pArray[_INDEX_chosen][k]*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
		                                +_DIS_z_eachatom.pArray[_INDEX_chosen][k]*_DIS_z_eachatom.pArray[_INDEX_chosen][k];
		_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]=0.0;
		for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
			for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
				for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
					_DIS_x=_PBL_X*x;
					_DIS_y=_PBL_Y*y;
					_DIS_z=_PBL_Z*z;
					_DIS2=_DIS2_eachatom.pArray[_INDEX_chosen][k]+2*_DIS_x*_DIS_x_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_y*_DIS_y_eachatom.pArray[_INDEX_chosen][k]
					                                      +2*_DIS_z*_DIS_z_eachatom.pArray[_INDEX_chosen][k]
					                                      +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
					_DIS6=TRIPLE(_DIS2);
					_DIS12=_DIS6*_DIS6;
					if(_DIS2<_R_cut_RR_eachatom.pArray[_INDEX_chosen][k]) {
						_ENER_LJ_eachatom.pArray[_INDEX_chosen][k]+=Energy_LJ(_EPSILON_eachatom.pArray[_INDEX_chosen][k],_SIGMA12_eachatom.pArray[_INDEX_chosen][k],_SIGMA6_eachatom.pArray[_INDEX_chosen][k],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[_INDEX_chosen][k];
					}
				}
			}
		}
	}
	//cout<<"]"<<endl;
	//cout<<" _ENER_eachatom["<<_INDEX_chosen<<"]="<<_ENER_eachatom[_INDEX_chosen]<<endl;
}
/////////////////////////
////////////////////////////
double cmc::calc_energy() {                       //the total energy;
	cout<<" calculating the total energy start (calc_energy)"<<endl;
	int i=0;
	int j=0;
	int k=0;
	int l=0;
	int x=0;
	int y=0;
	int z=0;

	_ENER_total=0.0;

	int Size_ATM=0;
	int Size_ATM_A=0;

	int index_atm=0;
	int index_atm_a=0; 

	for(i=0; i<_NUM_chains; i++) {// intra chain neighbor (hook<->chain) (L.J.<->others);
		Size_ATM=_SIZE_of_chn[i];
		if( _INDEX_chn_or_not_chain[i] ) {//if chain
			if( ( Size_ATM > 1 ) && _BF_flag ) {
				for(j=1; j<=Size_ATM-1; j++) {// must be '<', attention! 
				//index+1, cos' _XX _YY _ZZ defined to do this! 
					tempind_x=index_atm+j;
					tempind_y=index_atm+j+1;
					//cout<<tempind_x<<":"<<tempind_y<<endl;
					_DIS_x_eachatom.pArray[tempind_x][tempind_y]=_XX[tempind_x]-_XX[tempind_y];
					_DIS_y_eachatom.pArray[tempind_x][tempind_y]=_YY[tempind_x]-_YY[tempind_y];
					_DIS_z_eachatom.pArray[tempind_x][tempind_y]=_ZZ[tempind_x]-_ZZ[tempind_y];
					_DIS2_eachatom.pArray[tempind_x][tempind_y]=_DIS_x_eachatom.pArray[tempind_x][tempind_y]*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
				                                               +_DIS_y_eachatom.pArray[tempind_x][tempind_y]*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
				                                               +_DIS_z_eachatom.pArray[tempind_x][tempind_y]*_DIS_z_eachatom.pArray[tempind_x][tempind_y];
					_ENER_total+=Energy_BF(_PARA_K,_DIS2);
					//enumer++;
					//cout<<index_atm+j<<"<->"<<index_atm+j+1<<endl;
				}
			}
			if( ( Size_ATM > 2 ) && _AG_flag ) {
			}
			if( ( Size_ATM > 3 ) && _DH_flag ) {
			}
		} else {// if not chain 
			for(j=1; j<=Size_ATM-_NEIGHBOR_lj; j++) {//index+1, cos' _XX _YY _ZZ defined to do this! 
				for(k=1; k<=_NEIGHBOR_lj; k++) {
					tempind_x=index_atm+j;
					tempind_y=index_atm+j+k;
					_DIS_x_eachatom.pArray[tempind_x][tempind_y]=_XX[tempind_x]-_XX[tempind_y];
					_DIS_y_eachatom.pArray[tempind_x][tempind_y]=_YY[tempind_x]-_YY[tempind_y];
					_DIS_z_eachatom.pArray[tempind_x][tempind_y]=_ZZ[tempind_x]-_ZZ[tempind_y];
					_DIS2_eachatom.pArray[tempind_x][tempind_y]=_DIS_x_eachatom.pArray[tempind_x][tempind_y]*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
				                                        +_DIS_y_eachatom.pArray[tempind_x][tempind_y]*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
				                                        +_DIS_z_eachatom.pArray[tempind_x][tempind_y]*_DIS_z_eachatom.pArray[tempind_x][tempind_y];
					for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
						for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
							for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
								_DIS_x=_PBL_X*x;
								_DIS_y=_PBL_Y*y;
								_DIS_z=_PBL_Z*z;
								_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y]+2*_DIS_x*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_y*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_z*_DIS_z_eachatom.pArray[tempind_x][tempind_y]
								                                          +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
								_DIS6=TRIPLE(_DIS2);
								_DIS12=_DIS6*_DIS6;
								if(_DIS2<_R_cut_RR_eachatom.pArray[tempind_x][tempind_y]) {
									_ENER_total+=Energy_LJ(_EPSILON_eachatom.pArray[tempind_x][tempind_y],_SIGMA12_eachatom.pArray[tempind_x][tempind_y],_SIGMA6_eachatom.pArray[tempind_x][tempind_y],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[tempind_x][tempind_y];
								}
							}
						}
					}
					//enumer++;
				}
			}
		}
		//index_atm+=Size_ATM;
	}
	//cout<<" check point: "<<_ENER_total<<endl;
	index_atm=0;
	for(i=0; i<_NUM_chains; i++) {// intra chain, but not neighbor, ( LJ - all );
		Size_ATM=_SIZE_of_chn[i];
		for(j=1; j<=Size_ATM-_NEIGHBOR_lj-1; j++) {// must be '<', attention! 
            //index+1, cos' _XX _YY _ZZ defined to do this! 
			for(k=j+_NEIGHBOR_lj+1; k<=Size_ATM; k++) {// must be '<', attention! 
			//index+1, cos' _XX _YY _ZZ defined to do this! 
				tempind_x=index_atm+j;
				tempind_y=index_atm+k;
				_DIS_x_eachatom.pArray[tempind_x][tempind_y]=_XX[tempind_x]-_XX[tempind_y];
				_DIS_y_eachatom.pArray[tempind_x][tempind_y]=_YY[tempind_x]-_YY[tempind_y];
				_DIS_z_eachatom.pArray[tempind_x][tempind_y]=_ZZ[tempind_x]-_ZZ[tempind_y];
				_DIS2_eachatom.pArray[tempind_x][tempind_y]=_DIS_x_eachatom.pArray[tempind_x][tempind_y]*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
			                                        +_DIS_y_eachatom.pArray[tempind_x][tempind_y]*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
			                                        +_DIS_z_eachatom.pArray[tempind_x][tempind_y]*_DIS_z_eachatom.pArray[tempind_x][tempind_y];
				for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
					for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
						for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
							_DIS_x=_PBL_X*x;
							_DIS_y=_PBL_Y*y;
							_DIS_z=_PBL_Z*z;
							_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y]+2*_DIS_x*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
							                                          +2*_DIS_y*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
							                                          +2*_DIS_z*_DIS_z_eachatom.pArray[tempind_x][tempind_y]
							                                          +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
							_DIS6=TRIPLE(_DIS2); 
							_DIS12=_DIS6*_DIS6;
							if(_DIS2<_R_cut_RR_eachatom.pArray[tempind_x][tempind_y]) {
								_ENER_total+=Energy_LJ(_EPSILON_eachatom.pArray[tempind_x][tempind_y],_SIGMA12_eachatom.pArray[tempind_x][tempind_y],_SIGMA6_eachatom.pArray[tempind_x][tempind_y],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[tempind_x][tempind_y];
							}
						}
					}
				}
				//enumer++;
				//cout<<index_atm+j<<"<->"<<index_atm+k<<":";
				//printf("%20.15f\n", _ENER_total);
			}
		}
		index_atm+=Size_ATM;
	}

	index_atm=0;
	for(i=0; i<_NUM_chains-1; i++) {// inter chain, ( LJ - all );
		Size_ATM=_SIZE_of_chn[i];
		index_atm_a=0;
		for(j=0; j<i+1; j++) {
			index_atm_a+=_SIZE_of_chn[j];
		}
		for(j=i+1; j<_NUM_chains; j++) {// must be '<', attention! 
			Size_ATM_A=_SIZE_of_chn[j];
			for(k=1; k<=Size_ATM; k++) {// must be '<', attention! 
			//index+1, cos' _XX _YY _ZZ defined to do this!  	
				tempind_x=index_atm+k;
				for(l=1; l<=Size_ATM_A; l++) {//index+1, cos' _XX _YY _ZZ defined to do this!
					tempind_y=index_atm_a+l;
					_DIS_x_eachatom.pArray[tempind_x][tempind_y]=_XX[tempind_x]-_XX[tempind_y];
					_DIS_y_eachatom.pArray[tempind_x][tempind_y]=_YY[tempind_x]-_YY[tempind_y];
					_DIS_z_eachatom.pArray[tempind_x][tempind_y]=_ZZ[tempind_x]-_ZZ[tempind_y];
					_DIS2_eachatom.pArray[tempind_x][tempind_y]=_DIS_x_eachatom.pArray[tempind_x][tempind_y]*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
				                                        +_DIS_y_eachatom.pArray[tempind_x][tempind_y]*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
				                                        +_DIS_z_eachatom.pArray[tempind_x][tempind_y]*_DIS_z_eachatom.pArray[tempind_x][tempind_y];
					for(x=1-_PBC_Dim; x<_PBC_Dim; x++) {
						for(y=1-_PBC_Dim; y<_PBC_Dim; y++) {
							for(z=1-_PBC_Dim; z<_PBC_Dim; z++) {
								_DIS_x=_PBL_X*x;
								_DIS_y=_PBL_Y*y;
								_DIS_z=_PBL_Z*z;
								_DIS2=_DIS2_eachatom.pArray[tempind_x][tempind_y]+2*_DIS_x*_DIS_x_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_y*_DIS_y_eachatom.pArray[tempind_x][tempind_y]
								                                          +2*_DIS_z*_DIS_z_eachatom.pArray[tempind_x][tempind_y]
								                                          +_DIS_x*_DIS_x+_DIS_y*_DIS_y+_DIS_z*_DIS_z;
								_DIS6=TRIPLE(_DIS2);
								_DIS12=_DIS6*_DIS6;
								if(_DIS2<_R_cut_RR_eachatom.pArray[tempind_x][tempind_y]) {
									_ENER_total+=Energy_LJ(_EPSILON_eachatom.pArray[tempind_x][tempind_y],_SIGMA12_eachatom.pArray[tempind_x][tempind_y],_SIGMA6_eachatom.pArray[tempind_x][tempind_y],_DIS12,_DIS6)-_E_cut_RR_eachatom.pArray[tempind_x][tempind_y];
								}
							}
						}
					}				
					//enumer++;
				}
			}
			index_atm_a+=Size_ATM_A;
		}
		index_atm+=Size_ATM;
	}
	cout<<" calculating the total energy done! (calc_energy)"<<endl;
	return _ENER_total;
}
////////////////////////////////
void cmc::init_energy() {
    // bfflag_bak for calculating neighboring _DIS2_eachatom;
	double bfflag_bak=false;
	if(_BF_flag==false) {
		cout<<" temporarily change _BF_flag to: [true] (neighbor distance calculation)"<<endl;
		bfflag_bak=true;
		_BF_flag=true;
		_PARA_K=0.0;
	}

	cout<<" *--> energy initializing start .. "<<endl;	
	///////////////   code start  //////////////////
	//enumer=0;
	calc_energy();//necessary
	//cout<<" calc_energy() succ!"<<endl;
	//cout<<" --------------------- enumer = "<<enumer<<endl;
	//enumer=0;
	calc_energy_each_atom();
    double Total=get_energy();//necessary
	//cout<<" --------------------- enumer = "<<enumer<<endl;
	//cout<<" calc_energy_each_atom() succ!"<<endl;
	/* this is for the wall attraction */
	///////////////////// this is just for check at initial point!
	if( fabs(_ENER_total-Total) > 1e-8 ) {
		printf(" _ENER_total-total_ener=%20.15f\n", _ENER_total-Total);
		printf(" _ENER_total=%20.15f\n", _ENER_total);
		printf(" _ENER_get=%20.15f\n", Total);
		ErrorMSG("fabs(_ENER_total-total_ener) > 1e-8 ::cmc::init_energy");
		exit(LOGICERROR);
	}
    ///////////////   code start  //////////////////
    cout<<" *--> energy initializing done! "<<endl;

	if(bfflag_bak==true) {
		cout<<" change _BF_flag back to: [false] (neighbor distance calculation)"<<endl;
		_BF_flag=false;
	}
}
/////////////////////////////////
void cmc::calc_energy_each_atom() {
	for(int i=1; i<_SIZE_memo; i++) {
		make_choice(i);
		calc_energy_chosen_atom();
	}
	_Enery_initialized=true;
}
/////////////////////////////////
double cmc::get_energy() {
	double Ret_Ener=0.0;
	int i=0;
	int j=0;
	double Ret_Ener_LJ=0.0;
	for(i=1; i<_SIZE_memo; i++) {
		for(j=i+1; j<_SIZE_memo; j++) {
			Ret_Ener_LJ+=_ENER_LJ_eachatom.pArray[i][j];
		}
	}
	Ret_Ener+=Ret_Ener_LJ;
	////////////////////// this is just for check at initial point!
	return Ret_Ener;
}
void cmc::initialization() {
	load_parameters( "_parameters.pls" );
	load_conformation();
	write_parameters( "_parameters.pls" );
	memo_allocation();
	memo_setzero();
	memo_evaluation();
	init_epsilonsigma();//can not change the order;
    init_energy();
    chck_bond_len();//secure!!
    init_statistic();
}
////////////////////////////////////
///////////////////////////////////
/*void cmc::cout_parameters() {
	cout<<setiosflags(ios::left)<<endl;
	//cout<<setw(30)<<" ::_system_.natoms: "<<_system_.natoms<<endl;
	cout<<setw(30)<<" ::ee: "<<ee<<endl;
	cout<<setw(30)<<" ::_NUM_atoms: "<<_NUM_atoms<<endl;
	cout<<setw(30)<<" ::_totalener: "<<_ENER_total<<endl;
	cout<<setw(30)<<" ::_IFVERBOSE: "<<_IFVERBOSE<<endl;
	//cout<<setw(30)<<" ::_START_fromzero: "<<_START_fromzero<<endl;
	cout<<resetiosflags(ios::left)<<endl;
	//}
}*/
////////////////////////////////////
void cmc::check_energy() {
	double E_temp;
	double E_temp_2;
	E_temp=get_energy();
	E_temp_2=_ENER_total;
	calc_energy();
	printf(" T[%8.4f]: %20.15f-%20.15f(%20.15f)=%e\n", _TEMPERATURE, _ENER_total, E_temp, E_temp_2, _ENER_total-E_temp);//before fout_conf_PBC
	if( fabs(_ENER_total-E_temp)>1e-6) {
		cout<<" there's something wrong with energy, mpi termed!"<<endl;
		exit(LOGICERROR);
	}

}
////////////////////////////////////
void cmc::run_with_stepnumber(const int stepnumber) {

	cout<<endl<<" starting minimization run..... "<<endl<<endl;

	for(int i=0; i<stepnumber; i++) {
		make_mcmove( int(rand_seed(iseed_index)*_NUM_atoms)+1 );
	}
}
/////////////////
void cmc::run() {
	//int i=0;
	//int TEN_NUM_atoms=0;
	check_energy();

	cout<<endl<<" starting production run..... "<<endl<<endl;

	_MC_NUM_TOT=0;
	_MC_NUM_SUC=0;
	_MC_NUM_FIL=0;
	//int i=0;
	for(_I_eachstep=0, _I_totalnum=0; _I_totalnum < _RUNTIMES_totalnum; _I_eachstep++) {
		//recording(); //after every make_mcmove;
		if( _I_eachstep == _RUNTIMES_eachstep ) {
			if( (_I_totalnum+1)%_RUNTIMES_output == 0 ) {
				cout<<" _I_totalnum/_RUNTIMES_totalnum="<<setw(6)<<_I_totalnum+1
					<<"/"<<setw(6)<<setiosflags(ios::left)<<_RUNTIMES_totalnum
					<<resetiosflags(ios::left)<<endl;
			}
			_I_eachstep=0;
			_I_totalnum++;
			output_statistic();
			//trajectory_rec();			
		}
		make_mcmove( int(rand_seed(iseed_index)*_NUM_atoms)+1 );
		statistic();
	}
	check_energy();

	cout<<" succ rate : "<<double(_MC_NUM_SUC)/double(_MC_NUM_TOT)<<endl;
	cout<<" fail rate : "<<double(_MC_NUM_FIL)/double(_MC_NUM_TOT)<<endl;
		//fout_conformation_eachrep(true);
	//}
}
/////////////////////////////////
void cmc::init_statistic() {
	if(!_Enery_initialized) {
		cout<<" Energy not initialized, error: init_statistic() "<<endl;
		exit(LOGICERROR);
	}
	tcharn=new char[30];
	if(_FLAG_rh2) {
		_rh2_stream=new ofstream[_NUM_chains];
	}
	if(_FLAG_rg2) {
		_rg2_stream=new ofstream[_NUM_chains];
	}
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		stat_size=_SIZE_of_chn[stat_i];
		stat_head=_INDEX_CHN_HEAD[stat_i];
		stat_tail=_INDEX_CHN_TAIL[stat_i];
		if(_FLAG_rh2) {
			sprintf(tcharn,"rh2%03d%03d%05d.dat",(stat_i)+1,_NUM_chains,_NUM_atoms);
			_rh2_stream[stat_i].open(string(tcharn).c_str());
			cout<<" file: ["<<tcharn<<"] opened."<<endl;
			_rhead2[stat_i]=0.0;//_DIS2_eachatom.pArray[stat_head][stat_tail];... whatever...
		}
		if(_FLAG_rg2) {
			sprintf(tcharn,"rg2%03d%03d%05d.dat",(stat_i)+1,_NUM_chains,_NUM_atoms);
			_rg2_stream[stat_i].open(string(tcharn).c_str());
			cout<<" file: ["<<tcharn<<"] opened."<<endl;
			//initialization;
			_COM_x[stat_i]=0.0;
			_COM_y[stat_i]=0.0;
			_COM_z[stat_i]=0.0;
			for(stat_j=stat_head; stat_j<=stat_tail; stat_j++) {
				_COM_x[stat_i]+=_XX[stat_j];
				_COM_y[stat_i]+=_YY[stat_j];
				_COM_z[stat_i]+=_ZZ[stat_j];
			}
			//over
			_RG2_x[stat_i]=0.0;
			_RG2_y[stat_i]=0.0;
			_RG2_z[stat_i]=0.0;
		}
	}
}
////////////////
inline void cmc::statistic() { //every step!!
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		stat_size=_SIZE_of_chn[stat_i];
		stat_head=_INDEX_CHN_HEAD[stat_i];
		stat_tail=_INDEX_CHN_TAIL[stat_i];
		if(_FLAG_rh2) {
			_rhead2[stat_i]+=_DIS2_eachatom.pArray[stat_head][stat_tail];
		}
		if(_FLAG_rg2) {
			stat_com_x=_COM_x[stat_i]/stat_size;
			stat_com_y=_COM_y[stat_i]/stat_size;
			stat_com_z=_COM_z[stat_i]/stat_size;
			for(stat_j=stat_head; stat_j<=stat_tail; stat_j++) {
				tempnum=_XX[stat_j]-stat_com_x;
				_RG2_x[stat_i]+=tempnum*tempnum;
				tempnum=_YY[stat_j]-stat_com_y;
				_RG2_y[stat_i]+=tempnum*tempnum;
				tempnum=_ZZ[stat_j]-stat_com_z;
				_RG2_z[stat_i]+=tempnum*tempnum;
			}
		}
	}     
}
////////////////
inline void cmc::output_statistic() { //every _runtimes_eachstep
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		stat_size=_SIZE_of_chn[stat_i];
		stat_head=_INDEX_CHN_HEAD[stat_i];
		stat_tail=_INDEX_CHN_TAIL[stat_i];
		if(_FLAG_rh2) {
			_rh2_stream[stat_i]<<_rhead2[stat_i]/_RUNTIMES_eachstep<<endl;
			_rhead2[stat_i]=0.0;//_DIS2_eachatom.pArray[stat_head][stat_tail];
		}
		if(_FLAG_rg2) {
			_rg2_stream[stat_i]<<_RG2_x[stat_i]/_SIZE_of_chn[stat_i]/_RUNTIMES_eachstep<<" "
			                   <<_RG2_y[stat_i]/_SIZE_of_chn[stat_i]/_RUNTIMES_eachstep<<" "
			                   <<_RG2_z[stat_i]/_SIZE_of_chn[stat_i]/_RUNTIMES_eachstep<<endl;
			_RG2_x[stat_i]=0.0;
			_RG2_y[stat_i]=0.0;
			_RG2_z[stat_i]=0.0;
		}
	}       
}
////////////////
void cmc::close_statistic() { 
	char xtemp[30];
	sprintf(xtemp,"%d",_RUNTIMES_totalnum);
	for(stat_i=0; stat_i<_NUM_chains; stat_i++) {
		if(_FLAG_rh2) {
			sprintf(tcharn,"rh2%03d%03d%05d",(stat_i)+1,_NUM_chains,_NUM_atoms);
			system( (string("echo \"fn=\'") + tcharn + string("\';xlow=1;xhigh=")+string(xtemp)+string("\" > rh2.par")).c_str() );
			system("gnuplot < rh2.gpl");
			cout<<" file: ["<<tcharn<<".dat] closed, please check ["<<tcharn<<".eps]."<<endl;
			_rh2_stream[stat_i].close();
		}
		if(_FLAG_rg2) {
			sprintf(tcharn,"rg2%03d%03d%05d",(stat_i)+1,_NUM_chains,_NUM_atoms);
			system( (string("echo \"fn=\'") + tcharn + string("\';xlow=1;xhigh=")+string(xtemp)+string("\" > rg2.par")).c_str() );
			system("gnuplot < rg2.gpl");
			cout<<" file: ["<<tcharn<<".dat] closed, please check ["<<tcharn<<".eps]."<<endl;
			_rg2_stream[stat_i].close();
		}
	}   
	delete[] tcharn;
	if(_FLAG_rh2) {
		delete[] _rh2_stream;
	}
	if(_FLAG_rg2) {
		delete[] _rg2_stream;
	}
	_Statistic_over=true;
}
////////////////
void cmc::translate2box() {
	int i=0;
	for(i=1; i<_SIZE_memo; i++)
	{
		_XX[i]=Coor_PBC_X(_XX[i]);
		_YY[i]=Coor_PBC_Y(_YY[i]);
		_ZZ[i]=Coor_PBC_Z(_ZZ[i]);
	}
}
////////////////////////////////