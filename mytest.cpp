// CMSC 341 - Spring 2022 - Project 4
#include "dnadb.h"
#include <random>
#include <vector>
enum RANDOM { UNIFORMINT, UNIFORMREAL, NORMAL };
class Random {
public:
    Random(int min, int max, RANDOM type = UNIFORMINT, int mean = 50, int stdev = 20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL) {
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean, stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min, max);
        }
        else { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min, (double)max);
        }
    }
    void setSeed(int seedNum) {
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum() {
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if (m_type == NORMAL) {
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while (result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT) {
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum() {
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result * 100.0) / 100.0;
        return result;
    }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};
class Tester {
public:
    bool testInsertionNormalCase();
    bool testFindErrorCase();
    bool testFindFewNonCollidingKeys();
    bool testFindCollidingKeysWithoutRehash();
    bool testRemoveFewNonCollidingKeys();
    bool testRemoveCollidingKeysWithoutRehash();
    bool testRehashAfterInsertion();
    bool testRehashCompletionLoadFactor();
    bool testRehashAfterRemoval();
    bool testRehashCompletionDeleteRatio();
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main() {

    Tester tester;
    {
        // testing normal case for insert
        cout << "\n- Testing insert(...), normal case:\n";
        if (tester.testInsertionNormalCase() == true)
            cout << endl << "\t-> Normal case (Inserted DNA object at correct position) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Inserted DNA object at correct position) of DnaDb() failed!\n\n\n";
    }

    {
        // testing error case for find
        cout << "\n- Testing getDNA(...), error case:\n";
        if (tester.testFindErrorCase() == true)
            cout << endl << "\t-> Error case (Returned Empty DNA object) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Error case (Returned Empty DNA object) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for find with no colliding keys
        cout << "\n- Testing getDNA(...), normal case with no colliding keys:\n";
        if (tester.testFindFewNonCollidingKeys() == true)
            cout << endl << "\t-> Normal case (Found DNA objects) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Found DNA objects) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for find with colliding keys
        cout << "\n- Testing getDNA(...), normal case with colliding keys:\n";
        if (tester.testFindCollidingKeysWithoutRehash() == true)
            cout << endl << "\t-> Normal case (Found DNA objects) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Found DNA objects) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for remove with few non-colliding keys
        cout << "\n- Testing remove(...), normal case with few non-colliding keys:\n";
        if (tester.testRemoveFewNonCollidingKeys() == true)
            cout << endl << "\t-> Normal case (Removed DNA object) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Removed DNA object) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for remove with few colliding keys
        cout << "\n- Testing remove(...), normal case with few colliding keys:\n";
        if (tester.testRemoveCollidingKeysWithoutRehash() == true)
            cout << endl << "\t-> Normal case (Removed DNA object without triggering rehash) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Removed DNA object without triggering rehash) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for rehash(insert)
        cout << "\n- Testing insert(...), normal case for rehash:\n";
        if (tester.testRehashAfterInsertion() == true)
            cout << endl << "\t-> Normal case (Rehash triggered after decent number of objects inserted) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Rehash triggered after decent number of objects inserted) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for rehash(insert)
        cout << "\n- Testing insert(...), normal case for rehash due to load factor:\n";
        if (tester.testRehashCompletionLoadFactor() == true)
            cout << endl << "\t-> Normal case (Rehash completion due to load factor triggering) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Rehash completion due to load factor triggering) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for rehash(remove)
        cout << "\n- Testing remove(...), normal case triggering rehash:\n";
        if (tester.testRehashAfterRemoval() == true)
            cout << endl << "\t-> Normal case (Rehash triggered after decent number of objects removed) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Rehash triggered after decent number of objects removed) of DnaDb() failed!\n\n\n";
    }

    {
        // testing normal case for rehash(remove)
        cout << "\n- Testing remove(...), normal case triggering and completing rehash:\n";
        if (tester.testRehashCompletionDeleteRatio() == true)
            cout << endl << "\t-> Normal case (Rehash triggered and completed after decent number of objects removed) of DnaDb() passed!\n\n\n";
        else
            cout << endl << "\t-> Normal case (Rehash triggered and completed after decent number of objects removed) of DnaDb() failed!\n\n\n";
    }

    return 0;
}

bool Tester::testInsertionNormalCase() {

    bool check = true;
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    // checking whether all data are inserted
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++) {
        check = check && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }
    
    if (check) {
        return true;
    }
    return false;
}

bool Tester::testFindErrorCase() {

    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    DNA x = dnadb.getDNA("aaa", 1234);
    if (x.getSequence() == "") {
        return true;
    }
    return false;
}

bool Tester::testFindFewNonCollidingKeys() {

    bool check = true;
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 5; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    // checking whether all data are inserted
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++) {
        check = check && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }

    if (check) {
        return true;
    }
    return false;
}

bool Tester::testFindCollidingKeysWithoutRehash() {

    bool check = true;
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    // checking whether all data are inserted
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++) {
        check = check && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }

    if (check && !dnadb.m_isOldActivated) {
        return true;
    }
    return false;
}

bool Tester::testRemoveFewNonCollidingKeys() {

    DNA dnaObject;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 49; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        if (i == 17) {

            dnaObject = dataObj;
        }
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    // checking whether data is removed successfully
    if (dnadb.remove(dnaObject) && !dnadb.m_rehashTriggeredLambda) {
        return true;
    }
    return false;
}

bool Tester::testRemoveCollidingKeysWithoutRehash() {

    bool check = true;
    DNA dnaObject;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 5; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        if (i == 1) {

            dnaObject = dataObj;
        }
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    // checking whether data is removed successfully
    if (dnadb.remove(dnaObject) && !dnadb.m_rehashTriggeredLambda) {
        return true;
    }
    return false;
}

bool Tester::testRehashAfterInsertion() {

    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 98; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);

        if (dnadb.m_rehashTriggeredLambda) {
            return true;
        }
    }

    return false;
}

bool Tester::testRehashCompletionLoadFactor() {

    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 98; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    if (dnadb.m_hashCompleted) {
        return true;
    }
    return false;
}

bool Tester::testRehashAfterRemoval() {

    bool check = true;
    DNA list[80];
    DNA obj;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 98; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        if (i < 80) {
            list[i] = dataObj;
        }
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    for (int i = 0; i < 80; i++) {
        check = dnadb.remove(list[i]);
    }

    if (check && dnadb.m_rehashTriggeredRemove) {
        return true;
    }

    return false;
}

bool Tester::testRehashCompletionDeleteRatio() {

    bool check = true;
    DNA list[80];
    DNA obj;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i = 0; i < 98; i++) {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        // saving data for later use
        if (i < 80) {
            list[i] = dataObj;
        }
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }

    for (int i = 0; i < 80; i++) {
        check = dnadb.remove(list[i]);
    }

    if (check && dnadb.m_hashCompleted) {
        return true;
    }

    return false;
}

unsigned int hashCode(const string str) {
    unsigned int val = 0;
    const unsigned int thirtyThree = 33;  // magic number from textbook
    for (int i = 0; i < str.length(); i++)
        val = val * thirtyThree + str[i];
    return val;
}
string sequencer(int size, int seedNum) {
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0, 3);
    rndObject.setSeed(seedNum);
    for (int i = 0; i < size; i++) {
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}