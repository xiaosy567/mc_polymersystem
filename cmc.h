//cmc.h

#ifndef _CMC_H 
#define _CMC_H


#include "cmolecule.h"
#include "energy.h"
#include <sstream>
#define _OUTPUT_LEN  50
#define _CONVERGENCE 1e-12
#define _CONF_DIR {string("_conformations")}

////////////////////////////////////
template<typename T>
void writeparameter(ofstream &OFSName, const string ParaName, T Para) {
	OFSName<<setw(20)<<setiosflags(ios::left)<<ParaName.c_str()<<" "
		   <<setw(20)<<Para<<resetiosflags(ios::left)<<endl;
}
void readparameter(const string TempString, const string ParaName, bool &BPara);
void readparameter(const string TempString, const string ParaName, int &IPara);
void readparameter(const string TempString, const string ParaName, float &FPara);
void readparameter(const string TempString, const string ParaName, double &DPara);
void readparameter(const string TempString, const string ParaName, string &SPara);
/*void writeparameter(ofstream &OFSName, const string ParaName, bool BPrara);
void writeparameter(ofstream &OFSName, const string ParaName, int IPrara);
void writeparameter(ofstream &OFSName, const string ParaName, float FPrara);
void writeparameter(ofstream &OFSName, const string ParaName, double DPrara);
void writeparameter(ofstream &OFSName, const string ParaName, string SPrara);
*/
class cmc {
public:
	cmc();
	~cmc();
	void load_parameters(const string FILENAME_para);
	void write_parameters(const string FILENAME_para);
	//void make_dir();
	void init_conformation_a(const int NUM_atoms, const int NUM_chains, const double LEN_bond, const string FILENAME_conf);
	void initialization();
	///////////////////////////////////////////////////////////////////////////
	void load_conformation(); //'a' or 'f'	
	///////////////////////////////////////////////////////////////////////////
	void memo_allocation();
	bool _Memory_freed;
	void memo_setzero();    
	void memo_evaluation();
	inline void check_memo();
	void init_epsilonsigma();
	///////////////////////////////////////////////////////////////////////////
    double get_energy();                       //should = _ener_total
	double calc_energy();                      //the total energy; 
	void check_energy();
	void   init_energy();                      //doing before move;
	bool _Enery_initialized;
	void   calc_energy_each_atom();      
	inline void   calc_energy_chosen_atom();          //before & after each atom change; <1+<2
	void init_statistic();
	void close_statistic();
	bool _Statistic_over;
	inline void make_choice(const int INDEX_coor);
	//inline void prep_change(); //now embedded in make_choice;
	inline bool make_change();

	inline void  make_judge();
	inline void  make_judge_mc();

	inline void make_accept();
	inline void make_reject_all();	
	inline void make_reject_only_move();	

	inline int  make_mcmove(const int INDEX_coor);
	inline int  make_mcmove_mc(const int INDEX_coor);

	//inline void recording();
	//inline void trajectory_rec();
	void make_mapping(); // from _XX, _YY, _ZZ to the vector of atoms.
	bool chck_bond_len();
    void fout_conformation(const string FILENAME_conf, bool CHCK_BOND_LEN_OR_NOT);
	void memo_free();
	//////////////////////////////////////////////////////////////////////////
	CMolecule _system_;
	string _FILENAME_conf;
	bool _IFVERBOSE;
	
	int _NUM_atoms;
	int _NUM_chains;
	int _INDEX_chosen;
	double _ENER_total;
private:
	///////////////for test////////////////////////
	//int enumer;
	///////////////////////////////////////////////
	int _SIZE_memo;
	///////////////for bond fluctuation ///////////
	bool _BF_flag;
	bool _AG_flag;
	bool _DH_flag;
	double _PARA_K;
	double _BOND_delta;
	double _BOND_length;
	
	CMyArray<double> _EPSILON_eachatom;
	CMyArray<double> _SIGMA_eachatom;
	CMyArray<double> _SIGMA3_eachatom;
	CMyArray<double> _SIGMA6_eachatom;
	CMyArray<double> _SIGMA9_eachatom;
	CMyArray<double> _SIGMA12_eachatom;
	CMyArray<double> _E_cut_RR_eachatom;
	CMyArray<double> _R_cut_RR_eachatom;
	//CMyArray<int>    _IF_RR_eachatom;
	//int*             _IF_RR_calc_or_not;

	string _EPSILON_FN;
	string _SIGMA_FN;

	//for use in calc_ener_each_atom;
	double* cos_angle;
	double* sin_angle;
	int tempind;
	int tempind_x;
	int tempind_y;
	//for use in calc_ener_each_atom;

	//for judge
	//int tempindex_judge_bak;
	//int tempindex_judge;
	int tempnumer_ret;
	double TempEner;
	//for judge

	int     _NEIGHBOR_lj;

	//double  _epsilon_surf;
	//double  _sigma_surf;
	//double  _sigma_surf_1_6;
	
	//bool _ITERATION_OR_NOT;
	//bool _WHAM_EACH_OR_NOT;
	//bool _ENSEMBLE_AVERAGE;
	double  _CoorX1;
	double  _CoorX2;
	double  _CoorY1;
	double  _CoorY2;
	double  _CoorZ1;
	double  _CoorZ2;
	double  _PBL_X;
	double  _PBL_Y;
	double  _PBL_Z;
	int     _PBC_Dim;
	double  _DIS_x;
	double  _DIS_y;
	double  _DIS_z;
	////////////////////////////////////////
	template<class T>
	inline T Coor_PBC_X(T coor_x) {
		return coor_x>_CoorX1?( coor_x-_PBL_X*int((coor_x-_CoorX1)/_PBL_X) ) : ( coor_x-_PBL_X*int((coor_x-_CoorX2)/_PBL_X) );
	}
	template<class T>
	inline T Coor_PBC_Y(T coor_y) {
		return coor_y>_CoorY1?( coor_y-_PBL_Y*int((coor_y-_CoorY1)/_PBL_Y) ) : ( coor_y-_PBL_Y*int((coor_y-_CoorY2)/_PBL_Y) );
	}
	template<class T>
	inline T Coor_PBC_Z(T coor_z) {
		return coor_z>_CoorZ1?( coor_z-_PBL_Z*int((coor_z-_CoorZ1)/_PBL_Z) ) : ( coor_z-_PBL_Z*int((coor_z-_CoorZ2)/_PBL_Z) );
	}
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_X(T Delta_X) {
		T temp1=fabs( fabs(Delta_X)>_PBL_X?(Delta_X-_PBL_X*int(Delta_X/_PBL_X)):Delta_X );
		T temp2=temp1>(_PBL_X/2)?(_PBL_X-temp1):temp1;
		return temp2;
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_Y(T Delta_Y) {
		T temp1=fabs( fabs(Delta_Y)>_PBL_Y?(Delta_Y-_PBL_Y*int(Delta_Y/_PBL_Y)):Delta_Y );
		T temp2=temp1>(_PBL_Y/2)?(_PBL_Y-temp1):temp1;
		return temp2;
	}
	/////////////////////////////////////////
	////////////////////////////////////////
	template<class T>
	inline T DIS_PBC_Z(T Delta_Z) {
		T temp1=fabs( fabs(Delta_Z)>_PBL_Z?(Delta_Z-_PBL_Z*int(Delta_Z/_PBL_Z)):Delta_Z );
		T temp2=temp1>(_PBL_Z/2)?(_PBL_Z-temp1):temp1;
		return temp2;
	}
	/////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void translate2box();

	//int    _INDEX_TEMPERATURE;
	double _TEMPERATURE;
	double _TEMPERATURE_REP;

	//char* _OUTPUTDIR; 
	/* when rem, the following may need change;
	FILE* _OUTfile;
	char* _FNoutput;
	vector<string> _SAVElist;
	int _NUM_SAVElist;*/
	
	/*double* _PARA_beta;
	double* _PARA_alpha;
	double* _PARA_beta_new;
	double* _PARA_alpha_new;*/
////////////////////////////////
	CMyArray<double> _ENER_LJ_eachatom; //L-J potential;
	double* _ENER_LJ_chosenatom_backup;
	double* _ENER_BF_eachatom; //bond fluctuation;
	double  _ENER_BF_chosenatom_backup;
	/*CMyArray<double> _ENER_AG_eachatom; //bond fluctuation;
	double*  _ENER_AG_chosenatom_backup;
	CMyArray<double> _ENER_DH_eachatom; //dihedral;
	double*  _ENER_DH_chosenatom_backup;*/
	CMyArray<double> _DIS2_eachatom; // in one box
	double* _DIS2_chosenatom_backup; // in one box
	CMyArray<double> _DIS_x_eachatom; // in one box
	double* _DIS_x_chosenatom_backup; // in one box
	CMyArray<double> _DIS_y_eachatom; // in one box
	double* _DIS_y_chosenatom_backup; // in one box
	CMyArray<double> _DIS_z_eachatom; // in one box
	double* _DIS_z_chosenatom_backup; // in one box
	double* _COM_x;
	double* _COM_y;
	double* _COM_z;
	double* _RG2_x;
	double* _RG2_y;
	double* _RG2_z;

    char* tcharn;
	bool _FLAG_rg2;
	double* _rgyration2; //for each chain...
	ofstream* _rg2_stream;
	//////////////
	bool _FLAG_rh2;
	double* _rhead2;     //for each chain...
	ofstream* _rh2_stream;
	ofstream _ener_stream;
	inline void statistic();
	inline void output_statistic();
	int stat_i;
	int stat_j;
	int stat_k;
	int stat_size;
	int stat_head;
	int stat_tail; 
	double stat_com_x;
	double stat_com_y;
	double stat_com_z;
	double tempnum;
	int* _INDEX_CHN_HEAD; //for each chain...
	int* _INDEX_CHN_TAIL; //for each chain...
	//int _INDEX_chnhead;
	//int _INDEX_chntail;
	int* _INDEX_LNEIGHBOR;
	int* _INDEX_RNEIGHBOR;
	int _Index_lneighbor_ind;
	int _Index_rneighbor_ind;

	//double _DIS;
	double _DIS2;
	double _DIS3;
	double _DIS6;
	double _DIS9;
	double _DIS12;

	double  _ENER_delta;
	// the following 4 is for the chosen atom;
	int   _INDEX_chn_ind;           // for the chosen one;  a
	int   _SIZE_of_chosen_chn;      // for the chosen one;  b
	int   _INDEX_chn_or_not;        // for the chosen one;  c
	//int   _INDEX_atm_in_chn_ind;    // for the chosen one;  d
	int   _TYPE_atom_ind;           // for the chosen one;  e
	//void  memoset_fnode();
	int*  _INDEX_chn_atom;         // a
	int*  _SIZE_of_chn;            // b
	int*  _INDEX_chn_or_not_atom;  // c.1
	int*  _INDEX_chn_or_not_chain; // c.2
	int*  _INDEX_atm_in_chn;  // d
	int*  _TYPE_atom;         // e
	int*  _INDEX_RES_ATM;  // f.1
	//int*  _INDEX_RES_RES;  // f.2
	////////////////////////////////////////////////////////
	///////////////for rotation calculation ////////////////
	double rmin;
	double rmax;
	double ddelta;
	double r2min;
	double r2max;
	double ee;
	double ff;
	double dd;
	double dd_sqrt;
	double* ec;
	double* fc;
	double* dc;
	double* vecDC;
	double* coorRES;
	double factor;
	double* tc;
	double tl;
	double _PI_double;
	double _PI_half;
	double Omega;
	double THETA;
	double uc;
	double vc;
	double wc;
	double TempBX;
	double TempBY;
	double TempBZ;

	double sin_PHI;
	double cos_PHI;
	//double sin_THETA_B;
	//double cos_THETA_B;
	double Len_CB;
	double Len_CB_xy;
	double Len_CB_xy_2;
	double sin_Omega;
	double cos_Omega;
	double sin_THETA;
	double cos_THETA;
	/////////////// ends here //////////////////////////////
public:
	//void make_preparation();
	//void make_conformations();
	////////////////// mpisub //////////////////
	//void init_mpi(int pargc, char **pargv);              // init mpi or not;
	//void init_parameters();                            // init and bcast parameters;
	//void init_conformations();  //initialization;
	//void load_conformations();  //just load conformations;
	void run_with_stepnumber(const int stepnumber);
	void run();
	inline void   rand_update(const int procid);

private:

	double* _XX;
	double* _YY;
	double* _ZZ;
	/*double* _XX_rec;
	double* _YY_rec;
	double* _ZZ_rec;
	int _NUM_rec;*/

	int  _RUNTIMES_eachstep;         // 100000
    int  _RUNTIMES_totalnum;         // 10000 (10000 * 100000)
	int  _RUNTIMES_output;           //50

	//int  _RUNTIMES_iteration;
	//bool _START_fromzero;
	//int  _ITER;

	int _I_eachstep;
	int _I_totalnum;
	
 	////////////////  rand  ////////////////
	inline void   rand_init(const int procid);
	inline double rand_seed(int &iseed);
	//inline void   rand_update(const int procid);
	int *ir;
	int iy;
	int jj;

	double ran;

	int iseed_zero;
	int iseed_len1;
	int iseed_angle;
	int iseed_index;
	int iseed_rand;	

	double _MC_NUM_TOT; //total;
	double _MC_NUM_SUC; //succeed; ==1
	double _MC_NUM_FIL; //fail;    ==0|-1
	double _MC_NUM_STT; //statistic;
	//void fout_ratio();
	//void cout_parameters();

};

////////////////////////////////////
#define mm 714025
#define ia 1366
#define ic 150889
#define rm 1.0/mm
void cmc::rand_init(const int procid) {
	//different process has different seeds
	iseed_zero=(-2)*(procid+1)-1;
	iseed_len1=(-14)*(procid+1)-1;
	iseed_angle=(-6)*(procid+1)-1;
	iseed_index=(-10)*(procid+1)-1;
	iseed_rand=(-22)*(procid+1)-1;
}
void cmc::rand_update(const int procid) {
    iseed_zero=(-2)*(procid+1)-1;
	iseed_len1=(-14)*(procid+1)-1;
	iseed_angle=(-6)*(procid+1)-1;
	iseed_index=(-10)*(procid+1)-1;
	iseed_rand=(-22)*(procid+1)-1;
	MEMOSETZERO(ir, sizeof(int)*98);
	iy=0;
	jj=0;
}
double cmc::rand_seed(int &iseed) {
	if(iseed<0) {
		iseed=(int)fmod(double(ic-iseed), double(mm));
		for(int j=1; j<=97; j++) {
			iseed=(int)fmod(double(ia*iseed+ic), double(mm));
			ir[j]=iseed;
		}
		iseed=(int)fmod(double(ia*iseed+ic), double(mm));
		iy=iseed;
	}
	jj=1+(97*iy)/mm;
	iy=ir[jj];
	ran=iy*rm;
	iseed=(int)fmod(double(ia*iseed+ic), double(mm));
	ir[jj]=iseed;
	return ran;
}

/////////////////////////////////
#endif

