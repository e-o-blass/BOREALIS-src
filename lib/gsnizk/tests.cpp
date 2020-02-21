#define RUNS 100


/*
 * Copyright (c) 2016, Remi Bazin <bazin.remi@gmail.com>
 * See LICENSE for licensing details.
 */

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <chrono>
#include <sys/stat.h>

#include "gsnizk.h"
 
using namespace std;
using namespace gsnizk;

static int n_err = 0;
#define ASSERT(X) if (!(X)) (cerr << "Error: Assert of " << #X << " at line " \
    << __LINE__ << " failed!" << endl, ++n_err)

#define TRANSFER_TESTS 10
#define PAIRING_TESTS 10
#define PAIRING_COUNT_MAX 10
#define HASH_TESTS 10000

#define DATA_SIZE 2048

bool checkDataSize(int size, int line) {
    if (size > DATA_SIZE) {
        cerr << "Error: Please increase DATA_SIZE for the tests (size " << size
             << " required at line " << line << ")." << endl;
        pairings::terminate_pairings();
        return false;
    }
    return true;
}

#define CHECK_DATA_SIZE(s) if (!checkDataSize((s), __LINE__)) return

void testHash() {
    cout << "########## HASH TESTS ##########" << endl;
    srand(42);
    char *hash = new char[pairings::getHashLen()], *data = new char[256];
    ofstream out("hashes.test");
    for (int i = 0; i < HASH_TESTS; ++i) {
        int len = rand() % 257;
        for (int j = 0; j < len; ++j)
            data[j] = static_cast<char>(rand() & 0xFF);
        pairings::getHash(data, len, hash);
        out.write(hash, pairings::getHashLen());
    }
    out.close();
    /* Note: The file "hashes" is useful to check that all implementations
     * produce the same hashes. This is however limited to the getHash
     * function that uses SHA256 or SHA512. */
}



void testPairings() {
    cout << "########## PAIRING TESTS ##########" << endl;
    int len;
    char data[DATA_SIZE], *data2;

    /* -------------------- Hash prerequisite -------------------- */
    char *hash = new char[pairings::getHashLen()];
    pairings::getHash("hello", 5, hash);

    /* -------------------- Tests for Fp -------------------- */
    Fp v1, v2(0), v3(42), v4(1764);
    ASSERT(v1 == v2);
    v3 *= v3;
    ASSERT(v3 == v4);
    v2 += v3;
    ASSERT(v2 == v4);
    v2 -= v3;
    ASSERT(v1 == v2);
    v1 = Fp::getRand();
    ASSERT(v1 != v2); // Note: Just highly unlikely if the randomness is fine
    len = v1.getDataLen();
    cout << "Len for random Fp: " << len << endl;
    CHECK_DATA_SIZE(len);
    cout << "Testing transfers..." << endl;
    for (int i = 0; i < TRANSFER_TESTS; ++i) {
        v1 = Fp::getRand();
        v1.getData(data);
        v2 = Fp::getValue(data);
        ASSERT(v1 == v2);
    }
    v1 = Fp::getRand();
    v3 = (v1 / Fp(42)) * Fp(42);
    ASSERT(v1 == v3);
    ASSERT(Fp::fromHash("hello", 5) != Fp::fromHash("hi", 2));
    ASSERT(Fp::fromHash("hello", 5) == Fp::fromHash(hash));

    /* -------------------- Tests for G1 -------------------- */
    G1 g1 = G1::getRand(), g2, g3, g4;
    ASSERT(g1 != g2); // Note: Just highly unlikely if the randomness is fine
    len = g1.getDataLen();
    cout << "Len for random G1: " << len << endl;
    CHECK_DATA_SIZE(len);
    cout << "Testing transfers..." << endl;
    for (int i = 0; i < TRANSFER_TESTS; ++i) {
        g1 = G1::getRand();
        g1.getData(data);
        g2 = G1::getValue(data);
        ASSERT(g1 == g2);
    }
    len = g1.getDataLen(true);
    cout << "Len for random G1 compressed: " << len << endl;
    CHECK_DATA_SIZE(len);
    cout << "Testing transfers..." << endl;
    for (int i = 0; i < TRANSFER_TESTS; ++i) {
        g1 = G1::getRand();
        g1.getData(data, true);
        g3 = G1::getValue(data, true);
        ASSERT(g1 == g3);
    }
    g1 = G1::getRand();
    g3 = g1;
    g3 += g3;
    ASSERT((Fp(3) * g1) == (g1 + g3));
    ASSERT((-g1) == (Fp(-1) * g1));
    ASSERT((g1 - g1).isNull());
    ASSERT(G1::fromHash("hello", 5) != G1::fromHash("hi", 2));
    ASSERT(G1::fromHash("hello", 5) == G1::fromHash(hash));
    g1 = G1::getRand();
    v1 = Fp::getRand();
    g2 = v1 * g1;
    g3 = Fp::getUnit() * g1;
    ASSERT(g1 == g3);
    g1.precomputeForMult();
    g1.saveMultPrecomputations(data2);
    g3.loadMultPrecomputations(data2);
    ASSERT(g1 == g3);
    ASSERT(g2 == (v1 * g1));
    ASSERT(g2 == (v1 * g3));

    /* -------------------- Tests for G2 -------------------- */
    G2 h1 = G2::getRand(), h2, h3, h4;
    ASSERT(h1 != h2); // Note: Just highly unlikely if the randomness is fine
    len = h1.getDataLen();
    cout << "Len for random G2: " << len << endl;
    CHECK_DATA_SIZE(len);
    cout << "Testing transfers..." << endl;
    for (int i = 0; i < TRANSFER_TESTS; ++i) {
        h1 = G2::getRand();
        h1.getData(data);
        h2 = G2::getValue(data);
        ASSERT(h1 == h2);
    }
    len = h1.getDataLen(true);
    cout << "Len for random G2 compressed: " << len << endl;
    CHECK_DATA_SIZE(len);
    cout << "Testing transfers..." << endl;
    for (int i = 0; i < TRANSFER_TESTS; ++i) {
        h1 = G2::getRand();
        h1.getData(data, true);
        h2 = G2::getValue(data, true);
        ASSERT(h1 == h2);
    }
    h1 = G2::getRand();
    h3 = h1;
    h3 += h3;
    ASSERT((Fp(3) * h1) == (h1 + h3));
    ASSERT((-h1) == (Fp(-1) * h1));
    ASSERT((h1 - h1).isNull());
    ASSERT(G2::fromHash("hello", 5) != G2::fromHash("hi", 2));
    ASSERT(G2::fromHash("hello", 5) == G2::fromHash(hash));
    h1 = G2::getRand();
    v1 = Fp::getRand();
    h2 = v1 * h1;
    h3 = Fp::getUnit() * h1;
    ASSERT(h1 == h3);
    h1.precomputeForMult();
    h1.saveMultPrecomputations(data2);
    h3.loadMultPrecomputations(data2);
    ASSERT(h2 == (v1 * h1));
    ASSERT(h2 == (v1 * h3));

    /* -------------------- Tests for GT -------------------- */
    GT t1 = GT::getRand(), t2, t3, t4;
    ASSERT(t1 != t2); // Note: Just highly unlikely if the randomness is fine
    t2 = GT::getRand();
    ASSERT(t1 != t2); // Note: Just highly unlikely if the randomness is fine
    len = t1.getDataLen();
    cout << "Len for random GT: " << len << endl;
    CHECK_DATA_SIZE(len);
    cout << "Testing transfers..." << endl;
    for (int i = 0; i < TRANSFER_TESTS; ++i) {
        t1 = GT::getRand();
        t1.getData(data);
        t2 = GT::getValue(data);
        ASSERT(t1 == t2);
    }
    t1 = GT::getRand();
    t3 = t1;
    t3 *= t3;
    ASSERT((t1 ^ Fp(3)) == (t1 * t3));
    ASSERT((GT() / t1) == (t1 ^ Fp(-1)));
    ASSERT((t1 / t1).isUnit());
    t1 = GT::getRand();
    v1 = Fp::getRand();
    t2 = t1 ^ v1;
    t3 = t1 ^ Fp::getUnit();
    ASSERT(t1 == t3);
    t1.precomputeForPower();
    t1.savePowerPrecomputations(data2);
    t3.loadPowerPrecomputations(data2);
    ASSERT(t2 == (t1 ^ v1));
    ASSERT(t2 == (t3 ^ v1));

    /* -------------------- Pairing tests -------------------- */
    g1 = G1::getRand();
    h1 = G2::getRand();
    cout << "Testing simple pairings..." << endl;
    for (int i = 0; i < PAIRING_TESTS; ++i) {
        v1 = Fp::getRand();
        v2 = Fp::getRand();
        ASSERT(GT::pairing(v1 * g1, v2 * h1) ==
               (GT::pairing(g1, h1) ^ (v1 * v2)));
    }
    cout << "Testing multiple pairings..." << endl;
    for (int i = 0; i < PAIRING_TESTS; ++i) {
        std::vector< std::pair<G1,G2> > pairs;
        int n = rand() % PAIRING_COUNT_MAX;
        pairs.reserve(n);
        t1 = GT();
        while (n--) {
            g1 = G1::getRand();
            h1 = G2::getRand();
            t1 *= GT::pairing(g1, h1);
            pairs.push_back(std::pair<G1,G2>(g1, h1));
        }
        ASSERT(t1 == GT::pairing(pairs));
    }
    ASSERT((GT::pairing(g1, h1) * GT::pairing(-g1, h1)).isUnit());
    ASSERT((GT::pairing(g1, h1) * GT::pairing(g1, -h1)).isUnit());
    g1 = G1::getRand();
    h1 = G2::getRand();
    t1 = GT::pairing(g1, h1);
    h3 = Fp::getUnit() * h1;
    ASSERT(h1 == h3);
    h1.precomputeForPairing();
    h1.savePairingPrecomputations(data2);
    h3.loadPairingPrecomputations(data2);
    ASSERT(t1 == GT::pairing(g1, h1));
    ASSERT(t1 == GT::pairing(g1, h3));

    /* -------------------- iostream tests -------------------- */
    cout << "Testing iostream serialization..." << endl;
    v1 = Fp::getRand();
    v2 = Fp();
    g1 = G1::getRand();
    g2.clear();
    h1 = G2::getRand();
    h2.clear();
    t1 = GT::getRand();
    t2.clear();
    {
        ofstream out("pairings.test");
        out << v1 << v2;
        out << g1 << g2;
        out << h1 << h2;
        out << t1 << t2;
        out.close();
    }
    {
        ifstream in("pairings.test");
        in >> v3 >> v4;
        in >> g3 >> g4;
        in >> h3 >> h4;
        in >> t3 >> t4;
        in.close();
    }
    ASSERT(v1 == v3);
    ASSERT(v2 == v4);
    ASSERT(g1 == g3);
    ASSERT(g2 == g4);
    ASSERT(h1 == h3);
    ASSERT(h2 == h4);
    ASSERT(t1 == t3);
    ASSERT(t2 == t4);
}

//https://stackoverflow.com/questions/5840148/how-can-i-get-a-files-size-in-c/6039648#6039648
long GetFileSize(std::string filename)
{
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


long int testProof(NIZKProof &proof, ProofData &d, const CRS &crs, CRS *verif = 0) {
  long size = -1;
  //printf(" * Is proof zero knowledge? ");
  if (!proof.isZeroKnowledge()) {
   printf(" * !!!Proof is not zero-knowledge!!!\n");
}
    ASSERT(proof.verifySolution(d, crs));
    {
      //cout << " * Creating and writing proof..." << endl;
        ofstream out("proof.test");
	proof.writeProof(out, crs, d);
	out.close();
	size = GetFileSize("proof.test");
        
    }
    d.privFp.clear();
    d.privG1.clear();
    d.privG2.clear();
    {
      //cout << " * Reading and checking proof..." << endl;
        ifstream in("proof.test");
        if (verif) {
            ASSERT(proof.checkProof(in, *verif, d));
        } else {
            ASSERT(proof.checkProof(in, crs, d));
        }
        in.close();
    }
    if (!(proof.isZeroKnowledge() && crs.isSimulationReady()))
        return -1;
    {
      //cout << " * Creating and writing simulated proof..." << endl;
        ofstream out("proof-sim.test");
        proof.simulateProof(out, crs, d);
        out.close();
    }
    {
      //cout << " * Reading and checking simulated proof..." << endl;
        ifstream in("proof-sim.test");
        ASSERT(proof.checkProof(in, crs, d));
        in.close();
    }
    return size;
}


void OLDproveBit() {

  auto start = chrono::steady_clock::now();
  for (int i = 0;i<RUNS;i++) {
    CRS crs(false);
    ProofData d; 
    Fp beta(1);
    FpElement TILDEBETA = FpVar(0);
    d.privFp.push_back(beta);
    NIZKProof proof;
    proof.addEquation(TILDEBETA*TILDEBETA, FpUnit()*TILDEBETA);
    proof.endEquations();
    proof.dummyWrite(crs, d);
  }
  auto diff = chrono::steady_clock::now() - start;
  cout << "Bit: " << chrono::duration <double, milli> (diff).count()/RUNS << " ms" ;



  
  {
    CRS crs(false);
    ProofData d; 
    Fp beta(1);
    FpElement TILDEBETA = FpVar(0);
    d.privFp.push_back(beta);
    NIZKProof proof;
    proof.addEquation(TILDEBETA*TILDEBETA, FpUnit()*TILDEBETA);
    ASSERT(proof.endEquations());

    int size = testProof(proof, d, crs);
    cout << ", proof size: " <<size<<" Byte"<<endl;
  }
}


void proveDLOG() {
  std::vector <Fp> sks;
  std::vector <G1> pks;
  std::vector <CRS> crss;
  for (int i = 0;i<RUNS;i++) {
    CRS crs(false);
    crss.push_back(crs);

    Fp sk = Fp::getRand();
    sks.push_back(sk);
    G1 pk = sk * crs.getG1Base();
    pks.push_back(pk);
  }
  
  auto start = chrono::steady_clock::now();
  for (int i = 0;i<RUNS;i++) {
    
    ProofData d;
    FpElement TILDESK = FpVar(0);
    d.privFp.push_back(sks[i]);
    
    d.pubG1.push_back(pks[i]);
    G1Element PK = G1Const(0);

    G1 g1base = crss[i].getG1Base();
    d.pubG1.push_back(g1base);
    G1Element G1BASE = G1Const(1);
        
    NIZKProof proof;
    proof.addEquation(TILDESK * G1BASE, FpUnit()*PK);
    ASSERT(proof.endEquations());
    proof.dummyWrite(crss[i], d);
  }
  auto diff = chrono::steady_clock::now() - start;
  cout << "DLOG: " << chrono::duration <double, milli> (diff).count()/RUNS << " ms" ;

  {
CRS crs(false);

    ProofData d;
    Fp sk = Fp::getRand();
    FpElement TILDESK = FpVar(0);
    d.privFp.push_back(sk);
    G1 pk = sk * crs.getG1Base();
    d.pubG1.push_back(pk);
    G1Element PK = G1Const(0);

    G1 g1base = crs.getG1Base();
    d.pubG1.push_back(g1base);
    G1Element G1BASE = G1Const(1);
        
    NIZKProof proof;
    proof.addEquation(TILDESK * G1BASE, FpUnit()*PK);
    ASSERT(proof.endEquations());
    
    int size = testProof(proof, d, crs);
    cout << ", proof size: " << size << " Byte"<< endl;
  }
}

void proveDecryption() {
  std::vector <Fp> sks;
  std::vector <G1> pks;
  std::vector <G1> pk2s;
  std::vector <CRS> crss;
  for (int i = 0;i<RUNS;i++) {
    CRS crs(false);
    crss.push_back(crs);

    Fp sk = Fp::getRand();
    sks.push_back(sk);
    G1 pk = sk * crs.getG1Base();
    pks.push_back(pk);

    G1 pk2 = sk * pk;
    pk2s.push_back(pk2);
  }
  
  auto start = chrono::steady_clock::now();
  for (int i = 0;i<RUNS;i++) {
    
    ProofData d;
    FpElement TILDESK = FpVar(0);
    d.privFp.push_back(sks[i]);
    
    d.pubG1.push_back(pks[i]);
    G1Element PK = G1Const(0);

    G1 g1base = crss[i].getG1Base();
    d.pubG1.push_back(g1base);
    G1Element G1BASE = G1Const(1);

    d.pubG1.push_back(pk2s[i]);
    G1Element PK2 = G1Const(2);

    d.pubG1.push_back(pks[i]);
    G1Element G1BASE2 = G1Const(3);
    
    NIZKProof proof;
    proof.addEquation(TILDESK * G1BASE, FpUnit()*PK);
    proof.addEquation(TILDESK * G1BASE2, FpUnit()*PK2);
    ASSERT(proof.endEquations());
    proof.dummyWrite(crss[i], d);
  }

  auto diff = chrono::steady_clock::now() - start;
  cout << "Decryption: " << chrono::duration <double, milli> (diff).count()/RUNS << " ms" ;

  {
CRS crs(false);

    ProofData d;
    Fp sk = Fp::getRand();
    FpElement TILDESK = FpVar(0);
    d.privFp.push_back(sk);
    G1 pk = sk * crs.getG1Base();
    d.pubG1.push_back(pk);
    G1Element PK = G1Const(0);

    G1 g1base = crs.getG1Base();
    d.pubG1.push_back(g1base);
    G1Element G1BASE = G1Const(1);

    d.pubG1.push_back(pk);
    G1Element G1BASE2 = G1Const(2);

    G1 pk2 = sk * pk;
    d.pubG1.push_back(pk2);
    G1Element PK2 = G1Const(3);

    
    NIZKProof proof;
    proof.addEquation(TILDESK * G1BASE, FpUnit()*PK);
    proof.addEquation(TILDESK * G1BASE2, FpUnit()*PK2);

    ASSERT(proof.endEquations());
    
    int size = testProof(proof, d, crs);
    cout << ", proof size: " << size << " Byte"<< endl;
  }
}




void proveEnc() {
  std::vector <Fp> ras;
  std::vector <Fp> aones;
  std::vector <G1> lefts;
  std::vector <G1> rights;
  std::vector <CRS> crss;
  for (int i = 0;i<RUNS;i++) {
    CRS crs(false);
    crss.push_back(crs);

    Fp ra = Fp::getRand();
    ras.push_back(ra);
    
    G1 left = ra * crs.getG1Base();
    lefts.push_back(left);

    Fp aone = Fp::getRand();
    aones.push_back(aone);

    G1 right = ra * left + aone * crs.getG1Base();
    rights.push_back(right);
  }
  
  auto start = chrono::steady_clock::now();
  for (int i = 0;i<RUNS;i++) {
    
    ProofData d;
    FpElement TILDERA = FpVar(0);
    d.privFp.push_back(ras[i]);

    FpElement TILDEAONE = FpVar(1);
    d.privFp.push_back(aones[i]);
    
    d.pubG1.push_back(lefts[i]);
    G1Element LEFT = G1Const(0);

    G1 g1base = crss[i].getG1Base();
    d.pubG1.push_back(g1base);
    G1Element G1BASE = G1Const(1);

    d.pubG1.push_back(rights[i]);
    G1Element RIGHT = G1Const(2);
    
    NIZKProof proof;
    proof.addEquation(TILDERA * G1BASE, FpUnit()*LEFT);

    proof.addEquation(TILDERA * LEFT + TILDEAONE * G1BASE, FpUnit()* RIGHT);
    
    ASSERT(proof.endEquations());
    proof.dummyWrite(crss[i], d);
    //testProof(proof, d, crss[i]);    
  }
  auto diff = chrono::steady_clock::now() - start;
  cout << "Enc: " << chrono::duration <double, milli> (diff).count()/RUNS << " ms";

    {

    ProofData d;
    FpElement TILDERA = FpVar(0);
    d.privFp.push_back(ras[0]);

    FpElement TILDEAONE = FpVar(1);
    d.privFp.push_back(aones[0]);
    
    d.pubG1.push_back(lefts[0]);
    G1Element LEFT = G1Const(0);

    G1 g1base = crss[0].getG1Base();
    d.pubG1.push_back(g1base);
    G1Element G1BASE = G1Const(1);

    d.pubG1.push_back(rights[0]);
    G1Element RIGHT = G1Const(2);
    
    NIZKProof proof;
    proof.addEquation(TILDERA * G1BASE, FpUnit()*LEFT);

    proof.addEquation(TILDERA * LEFT + TILDEAONE * G1BASE, FpUnit()* RIGHT);
    
    ASSERT(proof.endEquations());    
    int  size = testProof(proof, d, crss[0]);
    cout << ", proof size: " << size <<" Byte"<<endl;
  }
  
}


void DGKleft(double *time, int *size) {


  CRS crs(false);
    
    G1 P = crs.getG1Base();
    
    Fp b(1);

    Fp a1(1);
    Fp a2(1);
    Fp two(2);

    std::vector <G1> clefts;
    std::vector <G1> bob_two_lefts;
    std::vector <G1> bob_one_lefts;

    
    for (int i=0;i<RUNS;i++) {
    Fp r = Fp::getRand();
    Fp sk = Fp::getRand();
    G1 pk = sk * P;
    G1 bob_one_left = r * P;
    bob_one_lefts.push_back(bob_one_left);
    
    G1 bob_one_right = r * pk + b * P;

    r = Fp::getRand();
    G1 bob_two_left = r * crs.getG1Base();
    bob_two_lefts.push_back(bob_two_left);
        
    G1 cleft = -bob_one_left + bob_two_left - two * a2 * bob_two_left;
    clefts.push_back(cleft);
    
    }

    auto start = chrono::steady_clock::now();
    for (int i=0;i<RUNS;i++) {
    NIZKProof proof;
    ProofData d; 

    int g1var = 0;
    int fpvar = 0;
    int g1const = 0;
    
    d.privG1.push_back(clefts[i]);
    G1Element TILDEC1LEFT = G1Var(g1var++);

    d.privFp.push_back(a2);
    FpElement TILDEA2 = FpVar(fpvar++);


    G1 gamma1 = two * bob_two_lefts[i];
    d.pubG1.push_back(gamma1);
    G1Element GAMMA1 = G1Const(g1const++);

    G1 t = - bob_one_lefts[i] + bob_two_lefts[i];
    d.pubG1.push_back(t);
    G1Element T = G1Const(g1const);
  
    
    proof.addEquation(FpUnit() * TILDEC1LEFT  + TILDEA2 * GAMMA1, FpUnit()*T);
    ASSERT(proof.endEquations());
    proof.dummyWrite(crs, d);
  }
    auto diff = chrono::steady_clock::now() - start;
    *time = chrono::duration <double, milli> (diff).count()/RUNS;

    
  {
    CRS crs(false);
    
    G1 P = crs.getG1Base();
    
    Fp b(1);

    Fp r = Fp::getRand();
    Fp sk = Fp::getRand();
    G1 pk = sk * P;
    G1 bob_one_left = r * P;
    G1 bob_one_right = r * pk + b * P;

    r = Fp::getRand();
    G1 bob_two_left = r * crs.getG1Base();
    G1 bob_two_right = r * pk + b * P;

    Fp a1(1);
    Fp a2(1);
    Fp two(2);
        
    G1 cleft = -bob_one_left + bob_two_left - two * a2 * bob_two_left;
    G1 cright = a1 * P - bob_one_right + P + a2*P + bob_two_right - two*a2*bob_two_right;

    NIZKProof proof;
    ProofData d; 

    int g1var = 0;
    int fpvar = 0;
    int g1const = 0;
    
    d.privG1.push_back(cleft);
    G1Element TILDEC1LEFT = G1Var(g1var++);

    d.privFp.push_back(a2);
    FpElement TILDEA2 = FpVar(fpvar++);


    G1 gamma1 = two * bob_two_left;
    d.pubG1.push_back(gamma1);
    G1Element GAMMA1 = G1Const(g1const++);

    G1 t = - bob_one_left + bob_two_left;
    d.pubG1.push_back(t);
    G1Element T = G1Const(g1const);
  
    
    proof.addEquation(FpUnit() * TILDEC1LEFT  + TILDEA2 * GAMMA1, FpUnit()*T);
    ASSERT(proof.endEquations());
    *size = testProof(proof, d, crs);
    }

}

void DGKright(double *time, int *size) {
    Fp a2(1);
    Fp a1(1);
    Fp two(2);

    CRS crs(false);  
    G1 P = crs.getG1Base();


    Fp sk = Fp::getRand();
    G1 pk = sk * P;
    Fp b(1);
       
    Fp r = Fp::getRand();
    G1 bob_one_right = r * pk + b * P;
    
    r = Fp::getRand();
    G1 bob_two_right = r * pk + b * P;

    G1 cright = a1 * P - bob_one_right + P + a2*P + bob_two_right - two*a2*bob_two_right;

    std::vector <G1> crights;
    std::vector <G1> bob_one_rights;
    std::vector <G1> bob_two_rights;

    
for (int i = 0;i<RUNS;i++) {
    Fp sk = Fp::getRand();
    G1 pk = sk * P;
    Fp b(1);
       
    Fp r = Fp::getRand();
    G1 bob_one_right = r * pk + b * P;
    bob_one_rights.push_back(bob_one_right);
    
    r = Fp::getRand();
    G1 bob_two_right = r * pk + b * P;
    bob_two_rights.push_back(bob_two_right);
    
    G1 cright = a1 * P - bob_one_right + P + a2*P + bob_two_right - two*a2*bob_two_right;
    crights.push_back(cright);
  
 }
    
    auto start = chrono::steady_clock::now();
    for (int i = 0;i<RUNS;i++) {
    NIZKProof proof;
    ProofData d; 

    int g1var = 0;
    int fpvar = 0;
    int g1const = 0;

      
    d.privFp.push_back(a2);
    FpElement TILDEA2 = FpVar(fpvar++);
    d.privG1.push_back(crights[i]);
    G1Element TILDECRIGHT = G1Var(g1var++); 

    d.privFp.push_back(a1);
    FpElement TILDEA1 = FpVar(fpvar++);

    d.pubG1.push_back(P);
    G1Element PP = G1Const(g1const++);
    
    d.pubG1.push_back(-bob_one_rights[i]);
    G1Element BOR = G1Const(g1const++);

    d.pubG1.push_back(bob_two_rights[i]);
    G1Element BTR = G1Const(g1const++);
    
    d.pubG1.push_back(-two*bob_two_rights[i]);
    G1Element TBTR = G1Const(g1const++);

    proof.addEquation(FpUnit()*TILDECRIGHT, TILDEA1 * PP + FpUnit()*BOR + FpUnit()*PP + TILDEA2 * PP + FpUnit()*BTR + TILDEA2 * TBTR);
    
    ASSERT(proof.endEquations());
    proof.dummyWrite(crs, d);
    }
    auto diff = chrono::steady_clock::now() - start;
    *time = chrono::duration <double, milli> (diff).count()/RUNS;

        
    {
    int g1var = 0;
    int fpvar = 0;
    int g1const = 0;

    NIZKProof proof;
    ProofData d; 

    d.privFp.push_back(a2);
    FpElement TILDEA2 = FpVar(fpvar++);
    d.privG1.push_back(cright);
    G1Element TILDECRIGHT = G1Var(g1var++); 

    d.privFp.push_back(a1);
    FpElement TILDEA1 = FpVar(fpvar++);

    d.pubG1.push_back(P);
    G1Element PP = G1Const(g1const++);
    
    d.pubG1.push_back(-bob_one_right);
    G1Element BOR = G1Const(g1const++);

    d.pubG1.push_back(bob_two_right);
    G1Element BTR = G1Const(g1const++);
    
    d.pubG1.push_back(-two*bob_two_right);
    G1Element TBTR = G1Const(g1const++);

    proof.addEquation(FpUnit()*TILDECRIGHT, TILDEA1 * PP + FpUnit()*BOR + FpUnit()*PP + TILDEA2 * PP + FpUnit()*BTR + TILDEA2 * TBTR);
    
    ASSERT(proof.endEquations());
    *size = testProof(proof, d, crs);
    }

}


void proveDGK() {
  double time = 0;
  int size = 0;
  double time2 = 0;
  int size2 = 0;

  
  DGKleft(&time, &size);
  DGKright(&time2, &size2);

  if ((time!=-1) &&(time2!=-1)) {
  time =  time+time2;
  cout << "DGK: "<< time<< " ms";
  }

  if ((size!=-1) &&(size2!=-1)) {
  size =  size+size2;
  cout << ", proof size: "<< size<< " Byte"<<endl;
  }

  
}


void proveBit() {
  CRS crs(false);
  G1 P = crs.getG1Base();

  Fp v(1);
  G1 base =  crs.getG1Base();

  std::vector <Fp> rs,sks;
  std::vector <G1> deltaones,deltazeros,czeros,cones,cprimezeros,cprimeones,pks;
  
  for (int i=0;i<RUNS;i++) {
    Fp r = Fp::getRand();
    rs.push_back(r);

    Fp sk = Fp::getRand();
    sks.push_back(sk);

    G1 czero = r * crs.getG1Base();
    czeros.push_back(czero);
    
    G1 cone = sk * czero + v * crs.getG1Base();
    cones.push_back(cone);

    G1 cprimezero = v * czero + r * crs.getG1Base();
    cprimezeros.push_back(cprimezero);

    G1 cprimeone = v * cone + r * sk * crs.getG1Base();
    cprimeones.push_back(cprimeone);
    
    pks.push_back(sk * crs.getG1Base());
    deltazeros.push_back(cprimezero - czero);
    deltaones.push_back(cprimeone - cone);

/*deltazeros.push_back(cprimezero - czero);
  deltaones.push_back(cprimeone - cone);*/
  }

  ProofData dprime;
  
  auto start = chrono::steady_clock::now();
  for (int i=0;i<RUNS;i++) {
  NIZKProof proof;
  ProofData d;

  FpElement TILDEV = FpVar(0);
  d.privFp.push_back(v);
  
  FpElement TILDER = FpVar(1);
  d.privFp.push_back(rs[i]);
 
  FpElement TILDESK = FpVar(2);
  d.privFp.push_back(sks[i]);

  d.pubG1.push_back(czeros[i]);
  G1Element CZERO = G1Const(0);

  d.pubG1.push_back(cones[i]);
  G1Element CONE = G1Const(1);
  
  d.pubG1.push_back(cprimezeros[i]);
  G1Element CPRIMEZERO = G1Const(2);

  d.pubG1.push_back(cprimeones[i]);
  G1Element CPRIMEONE = G1Const(3);

  d.pubG1.push_back(base);
  G1Element BASE = G1Const(4);

  d.pubG1.push_back(pks[i]);
  G1Element PK = G1Const(5);  

  d.pubG1.push_back(deltazeros[i]);
  G1Element DELTAZERO = G1Const(6);  

  d.pubG1.push_back(deltaones[i]);
  G1Element DELTAONE = G1Const(7);  

  
  proof.addEquation(TILDEV * CZERO + TILDER * BASE, FpUnit() * CPRIMEZERO);
  proof.addEquation(TILDEV * CONE + TILDER * PK, FpUnit() * CPRIMEONE);
  proof.addEquation(TILDESK * DELTAZERO, FpUnit() * DELTAONE);

  ASSERT(proof.endEquations());
  proof.dummyWrite(crs, d);
 
  }
   auto diff = chrono::steady_clock::now() - start;
   double time = chrono::duration <double, milli> (diff).count()/RUNS;
   cout << "Bit: " << time << " ms";

   {
   NIZKProof proof;
  ProofData d;

  FpElement TILDEV = FpVar(0);
  d.privFp.push_back(v);
  
  FpElement TILDER = FpVar(1);
  d.privFp.push_back(rs[0]);
 
  FpElement TILDESK = FpVar(2);
  d.privFp.push_back(sks[0]);

  d.pubG1.push_back(czeros[0]);
  G1Element CZERO = G1Const(0);

  d.pubG1.push_back(cones[0]);
  G1Element CONE = G1Const(1);
  
  d.pubG1.push_back(cprimezeros[0]);
  G1Element CPRIMEZERO = G1Const(2);

  d.pubG1.push_back(cprimeones[0]);
  G1Element CPRIMEONE = G1Const(3);

  d.pubG1.push_back(base);
  G1Element BASE = G1Const(4);

  d.pubG1.push_back(pks[0]);
  G1Element PK = G1Const(5);  

  d.pubG1.push_back(deltazeros[0]);
  G1Element DELTAZERO = G1Const(6);  

  d.pubG1.push_back(deltaones[0]);
  G1Element DELTAONE = G1Const(7);  

  
  proof.addEquation(TILDEV * CZERO + TILDER * BASE, FpUnit() * CPRIMEZERO);
  proof.addEquation(TILDEV * CONE + TILDER * PK, FpUnit() * CPRIMEONE);
  proof.addEquation(TILDESK * DELTAZERO, FpUnit() * DELTAONE);

  ASSERT(proof.endEquations());
  proof.dummyWrite(crs, d);
   printf(", proof size: %ld Byte\n",testProof(proof, d, crs));
}
  /*  Fp v(1);
  FpElement TILDEV = FpVar(0);
  d.privFp.push_back(v);

  Fp sk = Fp::getRand();
  FpElement TILDESK = FpVar(1);
  d.privFp.push_back(sk);

  Fp r = Fp::getRand();
  G1 czero = r * crs.getG1Base();
  d.pubG1.push_back(czero);
  G1Element CZERO = G1Const(0);

  G1 cone = sk * czero + v * crs.getG1Base();
  d.pubG1.push_back(cone);
  G1Element CONE = G1Const(1);
  
  G1 tildegammazero = czero - v * czero;
  d.privG1.push_back(tildegammazero);
  G1Element TILDEGAMMAZERO = G1Var(0);

  G1 tildegammaone = cone - v * cone;
  d.privG1.push_back(tildegammaone);
  G1Element TILDEGAMMAONE = G1Var(1);

  proof.addEquation(FpUnit() * TILDEGAMMAZERO + TILDEV * CZERO, FpUnit() * CZERO);
  proof.addEquation(FpUnit() * TILDEGAMMAONE + TILDEV * CONE, FpUnit() * CONE);
  proof.addEquation(TILDESK * TILDEGAMMAZERO, FpUnit() * TILDEGAMMAONE);

  ASSERT(proof.endEquations());
  proof.dummyWrite(crs, d);

  printf(", proof size: %ld Byte\n",testProof(proof, d, crs));
  */
}


void proveBlind() {
  CRS crs(false);
  G1 P = crs.getG1Base();

  std::vector <G1> crights;
  std::vector <G1> clefts;
    std::vector <G1> cprimerights;
      std::vector <G1> cprimelefts;
  std::vector <Fp> ras;

      
  for (int i=0;i<RUNS;i++) {
  Fp ra = Fp::getRand();
  ras.push_back(ra);
  G1 cleft = G1::getRand();
  clefts.push_back(cleft);
  G1 cright = G1::getRand();
  crights.push_back(cright);
  
  G1 cprimeleft = ra * cleft;
  cprimelefts.push_back(cprimeleft);
  G1 cprimeright = ra * cright;
  cprimerights.push_back(cprimeright);

  }

  auto start = chrono::steady_clock::now();
  for (int i=0;i<RUNS;i++) {
  NIZKProof proof;
  ProofData d;

  int g1var = 0;
  int fpvar = 0;

  d.privFp.push_back(ras[i]);
  FpElement TILDERA = FpVar(fpvar++);

  d.privG1.push_back(clefts[i]);
  G1Element TILDECLEFT = G1Var(g1var++);

  d.privG1.push_back(cprimelefts[i]);
  G1Element TILDECPRIMELEFT = G1Var(g1var++);


  d.privG1.push_back(crights[i]);
  G1Element TILDECRIGHT = G1Var(g1var++);

  d.privG1.push_back(cprimerights[i]);
  G1Element TILDECPRIMERIGHT = G1Var(g1var++);

  
  proof.addEquation(TILDERA * TILDECLEFT, FpUnit() * TILDECPRIMELEFT);
  proof.addEquation(TILDERA * TILDECRIGHT, FpUnit() * TILDECPRIMERIGHT);
  
  ASSERT(proof.endEquations());
  proof.dummyWrite(crs, d);
  }
   auto diff = chrono::steady_clock::now() - start;
   double time = chrono::duration <double, milli> (diff).count()/RUNS;
   cout << "Blind: " << time << " ms";
  
  {

  Fp ra = Fp::getRand();
  G1 cleft = G1::getRand();
  G1 cright = G1::getRand();
  
  G1 cprimeleft = ra * cleft;
  G1 cprimeright = ra * cright;

    
  NIZKProof proof;
  ProofData d;

  int g1var = 0;
  int fpvar = 0;

  d.privFp.push_back(ra);
  FpElement TILDERA = FpVar(fpvar++);

  d.privG1.push_back(cleft);
  G1Element TILDECLEFT = G1Var(g1var++);

  d.privG1.push_back(cprimeleft);
  G1Element TILDECPRIMELEFT = G1Var(g1var++);


  d.privG1.push_back(cright);
  G1Element TILDECRIGHT = G1Var(g1var++);

  d.privG1.push_back(cprimeright);
  G1Element TILDECPRIMERIGHT = G1Var(g1var++);

  
  proof.addEquation(TILDERA * TILDECLEFT, FpUnit() * TILDECPRIMELEFT);
  proof.addEquation(TILDERA * TILDECRIGHT, FpUnit() * TILDECPRIMERIGHT);
  
  ASSERT(proof.endEquations());
  printf(", proof size: %ld Byte\n",testProof(proof, d, crs));
  }
}

void proveShuffle() {

  CRS crs(false);
  Fp beta(1);

std::vector <G1> Clefts;
std::vector <G1> Clefttwos;
std::vector <G1> Crights;
 std::vector <G1> Crighttwos;
std::vector <G1> cprimetworights;
 std::vector <G1> cprimeonerights;
std::vector <G1> cprimetwolefts;
std::vector <G1> cprimeonelefts;
  for (int i=0;i<RUNS;i++)  {
  G1 cprimeoneleft = G1::getRand();
  cprimeonelefts.push_back(cprimeoneleft);
  
  G1 cprimetwoleft = G1::getRand();
  cprimetwolefts.push_back(cprimetwoleft);
  
  G1 cprimeoneright = G1::getRand();
  cprimeonerights.push_back(cprimeoneright);
  
  G1 cprimetworight = G1::getRand();
  cprimetworights.push_back(cprimetworight);

  G1 Cleft = beta * cprimeoneleft  + cprimetwoleft - beta * cprimetwoleft;
  Clefts.push_back(Cleft);
  
  G1 Cright = beta * cprimeoneright  + cprimetworight - beta * cprimetworight;
  Crights.push_back(Cright);

  G1 Clefttwo = beta * cprimetwoleft  + cprimeoneleft - beta * cprimeoneleft;
  Clefttwos.push_back(Clefttwo);
  
  G1 Crighttwo = beta * cprimetworight  + cprimeoneright - beta * cprimeoneright;
  Crighttwos.push_back(Crighttwo);
  }

  auto start = chrono::steady_clock::now();  
  for (int i=0;i<RUNS;i++) {
  
  NIZKProof proof;
  ProofData d;
  int g1const = 0;
  int fpvar = 0;
  int g1var = 0;
  
  d.pubG1.push_back(Clefts[i]);
  G1Element CLEFT = G1Const(g1const++);

  d.pubG1.push_back(Clefttwos[i]);
  G1Element CLEFTTWO = G1Const(g1const++);

  d.privFp.push_back(beta);
  FpElement TILDEBETA = FpVar(fpvar++);
  
  d.privFp.push_back(-beta);
  FpElement TILDEMINUSBETA = FpVar(fpvar++);

  d.privG1.push_back(cprimeonelefts[i]);
  G1Element TILDECPRIMEONELEFT = G1Var(g1var++);

  d.privG1.push_back(cprimetwolefts[i]);
  G1Element TILDECPRIMETWOLEFT = G1Var(g1var++);

  d.pubG1.push_back(Crights[i]);
  G1Element CRIGHT = G1Const(g1const++);

  d.pubG1.push_back(Crighttwos[i]);
  G1Element CRIGHTTWO = G1Const(g1const++);

  d.privG1.push_back(cprimeonerights[i]);
  G1Element TILDECPRIMEONERIGHT = G1Var(g1var++);

  d.privG1.push_back(cprimetworights[i]);
  G1Element TILDECPRIMETWORIGHT = G1Var(g1var++);

  
  proof.addEquation(TILDEBETA * TILDECPRIMEONELEFT  + FpUnit() * TILDECPRIMETWOLEFT + TILDEMINUSBETA * TILDECPRIMETWOLEFT, FpUnit() * CLEFT);

  proof.addEquation(TILDEBETA * TILDECPRIMEONERIGHT  + FpUnit() * TILDECPRIMETWORIGHT + TILDEMINUSBETA * TILDECPRIMETWORIGHT, FpUnit() * CRIGHT);


proof.addEquation(TILDEBETA * TILDECPRIMETWOLEFT  + FpUnit() * TILDECPRIMEONELEFT + TILDEMINUSBETA * TILDECPRIMEONELEFT, FpUnit() * CLEFTTWO);

  proof.addEquation(TILDEBETA * TILDECPRIMETWORIGHT  + FpUnit() * TILDECPRIMEONERIGHT + TILDEMINUSBETA * TILDECPRIMEONERIGHT, FpUnit() * CRIGHTTWO);
  
  ASSERT(proof.endEquations());
  proof.dummyWrite(crs, d);
}
   auto diff = chrono::steady_clock::now() - start;
   double time = chrono::duration <double, milli> (diff).count()/RUNS;
   cout << "Shuffle: " << time << " ms";
  
  {
  CRS crs(false);
  Fp beta(1);
  G1 cprimeoneleft = G1::getRand();
  G1 cprimetwoleft = G1::getRand();

  G1 cprimeoneright = G1::getRand();
  G1 cprimetworight = G1::getRand();

  G1 Cleft = beta * cprimeoneleft  + cprimetwoleft - beta * cprimetwoleft;
  G1 Cright = beta * cprimeoneright  + cprimetworight - beta * cprimetworight;


  G1 Clefttwo = beta * cprimetwoleft  + cprimeoneleft - beta * cprimeoneleft;
  G1 Crighttwo = beta * cprimetworight  + cprimeoneright - beta * cprimeoneright;

  
  NIZKProof proof;
  ProofData d;
  int g1const = 0;
  int fpvar = 0;
  int g1var = 0;
  
  d.pubG1.push_back(Cleft);
  G1Element CLEFT = G1Const(g1const++);

  d.pubG1.push_back(Clefttwo);
  G1Element CLEFTTWO = G1Const(g1const++);

  d.privFp.push_back(beta);
  FpElement TILDEBETA = FpVar(fpvar++);
  
  d.privFp.push_back(-beta);
  FpElement TILDEMINUSBETA = FpVar(fpvar++);

  d.privG1.push_back(cprimeoneleft);
  G1Element TILDECPRIMEONELEFT = G1Var(g1var++);

  d.privG1.push_back(cprimetwoleft);
  G1Element TILDECPRIMETWOLEFT = G1Var(g1var++);

  d.pubG1.push_back(Cright);
  G1Element CRIGHT = G1Const(g1const++);

  d.pubG1.push_back(Crighttwo);
  G1Element CRIGHTTWO = G1Const(g1const++);

  
  d.privG1.push_back(cprimeoneright);
  G1Element TILDECPRIMEONERIGHT = G1Var(g1var++);

  d.privG1.push_back(cprimetworight);
  G1Element TILDECPRIMETWORIGHT = G1Var(g1var++);

  
  proof.addEquation(TILDEBETA * TILDECPRIMEONELEFT  + FpUnit() * TILDECPRIMETWOLEFT + TILDEMINUSBETA * TILDECPRIMETWOLEFT, FpUnit() * CLEFT);

  proof.addEquation(TILDEBETA * TILDECPRIMEONERIGHT  + FpUnit() * TILDECPRIMETWORIGHT + TILDEMINUSBETA * TILDECPRIMETWORIGHT, FpUnit() * CRIGHT);


proof.addEquation(TILDEBETA * TILDECPRIMETWOLEFT  + FpUnit() * TILDECPRIMEONELEFT + TILDEMINUSBETA * TILDECPRIMEONELEFT, FpUnit() * CLEFTTWO);

  proof.addEquation(TILDEBETA * TILDECPRIMETWORIGHT  + FpUnit() * TILDECPRIMEONERIGHT + TILDEMINUSBETA * TILDECPRIMEONERIGHT, FpUnit() * CRIGHTTWO);

  
  
  ASSERT(proof.endEquations());
  printf(", proof size: %ld Byte\n",testProof(proof, d, crs));
  }
}

void testProofs() {
  //    cout << "########## PROOF TESTS ##########" << endl;

    CRS crs(false);
    CRS crsref(true), crspriv, crspub;
    crsref.makePublic();
    {
        ofstream out("crspriv.test");
        crspriv = crsref.genPrivate(out);
        out.close();
    }
    {
        ifstream in("crspriv.test");
        ASSERT(crsref.checkPrivate(in, crspriv));
        in.close();
    }
    remove("crspriv.test");
    crspub = crspriv;
    crspub.makePublic();
    /*
        {
        cout << "Instantiation 1: discrete log in G1" << endl;
        cout << " * Creating the equation system..." << endl;

        G1 a = G1::getRand();
        Fp k = Fp::getRand();
        G1 b = k * a;

	Fp myV(1);
	Fp myVV(2);
	G1 aa = a + a;
	
        NIZKProof proof;
	//original: proof.addEquation(FpVar(0) * G1Const(0), FpUnit() * G1Const(1));
	proof.addEquation(FpVar(0) * G1Const(0) + FpVar(1) * G1Const(1), FpVar(2) * G1Const(2));


	ASSERT(proof.endEquations());
        ProofData d;

			  
	d.privFp.push_back(k);
        d.pubG1.push_back(a);
	d.privFp.push_back(k);
        d.pubG1.push_back(a);
	d.privFp.push_back(myV);
        d.pubG1.push_back(b);
	
	//original
	//d.privFp.push_back(k);
        //.pubG1.push_back(a);
        //d.pubG1.push_back(b);
	
        testProof(proof, d, crs);
    }
    exit(1);
    {
        cout << "Instantiation 2: discrete log in G1 with private CRS" << endl;
        cout << " * Creating the equation system..." << endl;

        G1 a = G1::getRand();
        Fp k = Fp::getRand();
        G1 b = k * a;

        NIZKProof proof, proofcp;
        proof.addEquation(FpVar(0) * G1Const(a), FpUnit() * G1Const(b));

	ASSERT(proof.endEquations());

        ProofData d;
        d.privFp.push_back(k);
        

        cout << " * Writing and reading back the equation system..." << endl;
        {
            ofstream out("proof-model.test");
            out << proof;
            out.close();
        }
        {
            ifstream in("proof-model.test");
            in >> proofcp;
            in.close();
        }
        remove("proof-model.test");

        testProof(proofcp, d, crspriv, &crspub);
    }
    exit(0);*/
    {
      cout << "SCIB"<<endl<<RUNS <<" runs per operation" << endl;

        proveBit();
	proveDLOG();
	proveEnc();
	proveDGK();
	proveBlind();
	proveShuffle();
	proveDecryption();
	exit(0);
	
        ProofData d; 

	Fp sk = Fp::getRand();
	FpElement TILDESK = FpVar(0);
	d.privFp.push_back(sk);
	  
	G1 pk = sk * crs.getG1Base();
	d.pubG1.push_back(pk);
        G1Element PK = G1Const(0);

	G1 g1base = crs.getG1Base();
	d.pubG1.push_back(g1base);
        G1Element G1BASE = G1Const(1);

	
	G1 tildeg1base = crs.getG1Base();
	d.privG1.push_back(tildeg1base);
	G1Element TILDEG1BASE = G1Var(0);

	G1 tildeg1basetwo = crs.getG1Base();
	d.privG1.push_back(tildeg1basetwo);
	G1Element TILDEG1BASETWO = G1Var(1);
	
	Fp beta(1);
	FpElement TILDEBETA = FpVar(1);
	d.privFp.push_back(beta);
	
        NIZKProof proof;
	//secret Fp, public G1
	proof.addEquation(TILDESK * G1BASE, FpUnit()*PK);

	//public Fp, secret G1
	proof.addEquation(FpUnit() * TILDEG1BASE, FpUnit()*G1BASE);

	//quadratic in Fp
	proof.addEquation(TILDEBETA*TILDEBETA, FpUnit()*TILDEBETA);

	//mixing secrets in Fp and G1
	proof.addEquation(TILDESK * G1BASE + FpUnit() * TILDEG1BASE, FpUnit()*PK + FpUnit()*G1BASE);

	//multiple secrets in G1
	proof.addEquation(FpUnit() * TILDEG1BASE + FpUnit() * TILDEG1BASETWO, FpUnit()*G1BASE+FpUnit()*G1BASE);
	// Gamma != 0
	proof.addEquation(TILDESK * TILDEG1BASE, FpUnit()*PK);

	// Gamma != 0 and RHS contains variables
	proof.addEquation(TILDESK * TILDEG1BASE + TILDESK * TILDEG1BASE, FpUnit()*PK + TILDESK * TILDEG1BASE);
	
        ASSERT(proof.endEquations());

        testProof(proof, d, crs);
    }
    exit(0);
    
    {
        cout << "Instantiation 4: user tokens (2)" << endl;
        cout << "  (see https://eprint.iacr.org/2016/416)" << endl;
        cout << " * Creating the equation system..." << endl;

        ProofData d;

        /* Certificate Authority's credentials */
        Fp sk_A = Fp::getRand();
        G2 pk_A = sk_A * crs.getG2Base();
        d.pubG2.push_back(pk_A);
        G2Element _pk_A = G2Const(0);

        /* User's credentials */
        Fp sk_C = Fp::getRand();
        G1 pk_C = sk_C * crs.getG1Base();
        d.privFp.push_back(sk_C);
        FpElement _sk_C = FpVar(0);
        d.privG1.push_back(pk_C);
        G1Element _pk_C = G1Var(0);

        /* User's certificate */
        G1 cert = sk_A * pk_C;
        d.privG1.push_back(cert);
        G1Element _cert = G1Var(1);

        /* Hash of the one-time public key
         * (in G1 instead of G2 for efficiency) */
        G1 HK = G1::getRand();
        d.pubG1.push_back(HK);
        G1Element _HK = G1Const(0);

        /* Signature of the one-time public key */
        G1 sign = sk_C * HK;
        d.pubG1.push_back(sign);
        G1Element _sign = G1Const(1);

        /* Service Provider ID */
        G1 v_SP = G1::getRand();
        d.pubG1.push_back(v_SP);
        G1Element _v_SP = G1Const(2);

        /* Linkability value */
        G1 value = sk_C * v_SP;
        d.pubG1.push_back(value);
        G1Element _value = G1Const(3);

        NIZKProof proof;
        proof.addEquation(FpUnit() * _pk_C, _sk_C * G1Base());
        proof.addEquation(e(_cert, G2Base()), e(_pk_C, _pk_A));
        proof.addEquation(FpUnit() * _sign, _sk_C * _HK);
        proof.addEquation(FpUnit() * _value, _sk_C * _v_SP);
        ASSERT(proof.endEquations());

        testProof(proof, d, crs);
    }
    {
        cout << "Instantiation 5: Big equation" << endl;
        ProofData d;

        Fp k = Fp::getRand(), l = Fp::getRand();
        d.privFp.push_back(k);
        FpElement _k = FpVar(0);
        G1 v = (k * l) * crs.getG1Base();

        NIZKProof proof;
        proof.addEquation(e(G1Base(), (_k * FpConst(l)) * G2Base()),
                          e(G1Const(v), G2Base()));
        ASSERT(proof.endEquations());

        testProof(proof, d, crs);
    }
    {
        cout << "Instantiation 6: Extractable proof" << endl;
        ProofData d;

        CRS crs_extract(true);

        Fp k = Fp::getRand();
        G1 kg1 = k * crs_extract.getG1Base();
        G2 kg2 = k * crs_extract.getG2Base();

        G1Element _kg1 = G1Var(0);
        d.privG1.push_back(kg1);

        NIZKProof proof;
        proof.addEquation(e(_kg1, G2Base()), e(G1Base(), G2Const(kg2)));
        ASSERT(proof.endEquations());

        testProof(proof, d, crs_extract);

        cout << " * Extracting private value..." << endl;
        {
            ifstream in("proof.test");
            B1 c_kg1;
            in >> c_kg1;
            in.close();
            G1 recovered_kg1 = c_kg1.extract(crs_extract);
            ASSERT(recovered_kg1 == kg1);
        }
    }
    remove("proof.test");
    remove("proof-sim.test");
}

void testLibrary() {
  //testHash();
  //testPairings();
    testProofs();
    if (n_err) {
        cout << "Done; " << n_err << " error(s) have occured!" << endl;
    } else {
        cout << "Done; no errors have occured." << endl;
    }
}
