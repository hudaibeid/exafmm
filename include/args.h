#ifndef args_h
#define args_h
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#ifndef _SX
#include <getopt.h>
#endif

namespace exafmm {
#ifndef _SX
static struct option long_options[] = {
	{"ncrit",        required_argument, 0, 'c'},
	{"distribution", required_argument, 0, 'd'},
	{"dual",         no_argument,       0, 'D'},
	{"graft",        no_argument,       0, 'g'},
	{"getMatrix",    no_argument,       0, 'G'},
	{"help",         no_argument,       0, 'h'},
	{"images",       required_argument, 0, 'i'},
	{"IneJ",         no_argument,       0, 'j'},
	{"mutual",       no_argument,       0, 'm'},
	{"numBodies",    required_argument, 0, 'n'},
	{"useRopt",      no_argument,       0, 'o'},
	{"repeat",       required_argument, 0, 'r'},
	{"nspawn",       required_argument, 0, 's'},
	{"theta",        required_argument, 0, 't'},
	{"threads",      required_argument, 0, 'T'},
	{"verbose",      no_argument,       0, 'v'},
	{"write",        no_argument,       0, 'w'},
	{"useRmax",      no_argument,       0, 'x'},
	{"balance",      required_argument, 0, 'b'},
  {"weight",       required_argument, 0, 'W'},
  {"grain",        required_argument, 0, 'a'},
	{0, 0, 0, 0}
};
#endif

class Args {
public:
	int ncrit;
	const char * distribution;
	int dual;
	int graft;
	int getMatrix;
	int images;
	int IneJ;
	int mutual;
	int numBodies;
	int useRopt;
	int repeat;
	int nspawn;
	double theta;
	int threads;
	int verbose;
	int write;
	int useRmax;
	int balance;
	double weight;
	uint16_t grain;

private:
	void usage(char * name) {
		fprintf(stderr,
		        "Usage: %s [options]\n"
		        "Long option (short option)     : Description (Default value)\n"
		        " --ncrit (-c)                  : Number of bodies per leaf cell (%d)\n"
		        " --distribution (-d) [l/c/s/p] : lattice, cube, sphere, octant, plummer (%s)\n"
		        " --dual (-D)                   : Use dual tree traversal (%d)\n"
		        " --graft (-g)                  : Graft remote trees to global tree (%d)\n"
		        " --getMatrix (-G)              : Write G matrix to file (%d)\n"
		        " --help (-h)                   : Show this help document\n"
		        " --images (-i)                 : Number of periodic image levels (%d)\n"
		        " --IneJ (-j)                   : Use different sources & targets (%d)\n"
		        " --mutual (-m)                 : Use mutual interaction (%d)\n"
		        " --numBodies (-n)              : Number of bodies (%d)\n"
		        " --useRopt (-o)                : Use error optimized theta for MAC (%d)\n"
		        " --repeat (-r)                 : Number of iteration loops (%d)\n"
		        " --nspawn (-s)                 : Threshold for stopping task creation during recursion (%d)\n"
		        " --theta (-t)                  : Multipole acceptance criterion (%f)\n"
		        " --threads (-T)                : Number of threads (%d)\n"
		        " --verbose (-v)                : Print information to screen (%d)\n"
		        " --write (-w)                  : Write timings to file (%d)\n"
		        " --useRmax (-x)                : Use maximum distance for MAC (%d)\n"
		        " --balance (-b)                : Loadbalance input (%d)\n"
		        " --weight (-W)                 : remote weight (%f)\n"
		        " --grain (-a)                  : grain size for remote communication (%d)\n" ,
		        name,
		        ncrit,
		        distribution,
		        dual,
		        graft,
		        getMatrix,
		        images,
		        IneJ,
		        mutual,
		        numBodies,
		        useRopt,
		        repeat,
		        nspawn,
		        theta,
		        threads,
		        verbose,
		        write,
		        useRmax,
		        balance,
		        weight,
		        grain );
	}

	const char * parse(const char * arg) {
		switch (arg[0]) {
		case 'c':
			return "cube";
		case 'l':
			return "lattice";
		case 'o':
			return "octant";
		case 'p':
			return "plummer";
		case 's':		
			return "sphere";			
		case 'b':
			return "boundary";
		default:
			fprintf(stderr, "invalid distribution %s\n", arg);
			exit(0);
		}
		return "";
	}

public:
	Args(int argc = 0, char ** argv = NULL) :
		ncrit(64),
		distribution("cube"),
		dual(0),
		graft(0),
		getMatrix(0),
		images(0),
		IneJ(0),
		mutual(0),
		numBodies(1000000),
		useRopt(0),
		repeat(1),
		nspawn(5000),
		theta(.4),
		threads(16),
		verbose(0),
		write(0),
		useRmax(0),
		balance(0),
		weight(1.0),
		grain(1)
	{
		while (1) {
#if _SX
#warning SX does not have getopt_long
			int c = getopt(argc, argv, "c:d:DgGhi:jmn:or:s:t:T:b:W:a:vwx");
#else
			int option_index;
			int c = getopt_long(argc, argv, "c:d:DgGhi:jmn:or:s:t:T:b:W:a:vwx", long_options, &option_index);
#endif
			if (c == -1) break;
			switch (c) {
			case 'c':
				ncrit = atoi(optarg);
				break;
			case 'd':
				distribution = parse(optarg);
				break;
			case 'D':
				dual = 1;
				break;
			case 'g':
				graft = 1;
				break;
			case 'G':
				getMatrix = 1;
				break;
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'i':
				images = atoi(optarg);
				break;
			case 'j':
				IneJ = 1;
				break;
			case 'm':
				mutual = 1;
				break;
			case 'n':
				numBodies = atoi(optarg);
				break;
			case 'o':
				useRopt = 1;
				break;
			case 'r':
				repeat = atoi(optarg);
				break;
			case 's':
				nspawn = atoi(optarg);
				break;
			case 't':
				theta = atof(optarg);
				break;
			case 'T':
				threads = atoi(optarg);
				break;
			case 'v':
				verbose = 1;
				break;
			case 'w':
				write = 1;
				break;
			case 'x':
				useRmax = 1;
				break;
			case 'b':
				balance = atoi(optarg);
				break;
			case 'W':
				weight = atof(optarg);
				break;
			case 'a':
				grain = atoi(optarg);
				break;
			default:
				usage(argv[0]);
				exit(0);
			}
		}
	}
	template <typename OutputStream>
	void print(int stringLength, int PP, OutputStream&& stream) {
		if (verbose) {
			stream << std::setw(stringLength) << std::fixed << std::left
			          << "ncrit" << " : " << ncrit << std::endl
			          << std::setw(stringLength)
			          << "distribution" << " : " << distribution << std::endl
			          << std::setw(stringLength)
			          << "dual" << " : " << dual << std::endl
			          << std::setw(stringLength)
			          << "graft" << " : " << graft << std::endl
			          << std::setw(stringLength)
			          << "getMatrix" << " : " << getMatrix << std::endl
			          << std::setw(stringLength)
			          << "images" << " : " << images << std::endl
			          << std::setw(stringLength)
			          << "IneJ" << " : " << IneJ << std::endl
			          << std::setw(stringLength)
			          << "mutual" << " : " << mutual << std::endl
			          << std::setw(stringLength)
			          << "numBodies" << " : " << numBodies << std::endl
			          << std::setw(stringLength)
			          << "useRopt" << " : " << useRopt << std::endl
			          << std::setw(stringLength)
			          << "P" << " : " << PP << std::endl
			          << std::setw(stringLength)
			          << "repeat" << " : " << repeat << std::endl
			          << std::setw(stringLength)
			          << "nspawn" << " : " << nspawn << std::endl
			          << std::setw(stringLength)
			          << "theta" << " : " << theta << std::endl
			          << std::setw(stringLength)
			          << "threads" << " : " << threads << std::endl
			          << std::setw(stringLength)
			          << "verbose" << " : " << verbose << std::endl
			          << std::setw(stringLength)
			          << "write" << " : " << write << std::endl
			          << std::setw(stringLength)
			          << "useRmax" << " : " << useRmax << std::endl
    						<< std::setw(stringLength)                      
			          << "balance" << " : " << balance << std::endl   
    						<< std::setw(stringLength)
    						<< "remote weight" << " : " << weight << std::endl
						    << std::setw(stringLength)
						    << "comm grain" << " : " << int(grain) << std::endl;
		}
	}
};
}
#endif
