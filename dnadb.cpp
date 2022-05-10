// CMSC 341 - Spring 2022 - Project 4
#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash) {

    m_hash = hash;
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    
    if (size) {
        
        if (size > MAXPRIME) {
            m_currentCap = MAXPRIME;
        }
        else if (size < MINPRIME) {
            m_currentCap = MINPRIME;
        }
        else {

            // 2 cases: first, if it isnt a prime number, second, if it is.
            // if it isnt prime, find a prime which is greater than the immediate number
            if (isPrime(size)) {
                m_currentCap = size;
            }
            else {
                m_currentCap = findNextPrime(size);
            }
        }
    }

    m_currentTable = new DNA[m_currentCap];
    //m_oldTable = new DNA[m_currentCap];

    for (unsigned int i = 0; i < m_currentCap; i++) {
        m_currentTable[i] = EMPTY;
        //m_oldTable[i] = EMPTY;
    }
}

DnaDb::~DnaDb() {

}

bool DnaDb::insert(DNA dna) {

    int index = m_hash(dna.getSequence()) % m_currentCap;
    bool isInserted = false;

    if (m_currentTable[index] == EMPTY && isValidLocation(dna)) {

        m_currentTable[index] = dna;
        ++m_currentSize; // SAVING THE COUNT OF DNA DATA INSERTED
        isInserted = true;
    }
    else {

        // COLLISION OCCURED, quadratic probe to the next free location
        int count = 0;
        bool isAdded = false;
        int newIndexHash = m_hash(dna.getSequence()) % m_currentCap; 

        while (!isAdded) {

            index = (newIndexHash + (count * count)) % m_currentCap;
            if (m_currentTable[index] == EMPTY && isValidLocation(dna)) {

                m_currentTable[index] = dna;
                ++m_currentSize;
                isInserted = true;
                isAdded = true;
            }
            ++count;
        }
    }

    if (lambda() > 0.5) {

        //rehash the entire table
        m_isOldActivated = true;
        m_oldCap = m_currentCap;
        m_oldSize = m_currentSize;
        m_oldNumDeleted = m_currNumDeleted;

        m_oldTable = new DNA[m_oldCap];

        for (unsigned int i = 0; i < m_currentCap; i++) {
            
            m_oldTable[i] = m_currentTable[i];
        }

        delete[] m_currentTable;
        m_currentCap = (m_currentSize - m_currNumDeleted) * 4; //calculate the size
        m_currentSize = 0; // CHECK THIS JUST INCASE
        m_currentTable = new DNA[m_currentCap];

        
        int positionNewTable = 0; // INDEX FOR THE NEW TABLE
        switch (m_quarterChecker)
        {
        case 1:
            m_location = 0; // INITIALIZING THE LOCATION VARIABLE BECAUSE IS THE FIRST QUARTER OF INITIALIZATION
            m_quarterSize = int(ceil(m_oldSize / 4)); // STORES THE QUARTER SIZE BASED ON THE TABLE SIZE

            for (; m_location < m_quarterSize; m_location++) {

                if (m_oldTable[m_location].getSequence() != "" && m_oldTable[m_location].getSequence() != DELETEDKEY) {

                    positionNewTable = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;

                    // CHECK THIS JUST INCASE OF AN ERROR CAUSE I DIDNT INITIALIZE IT WITH THE EMPTY DNA OBJECT
                    if (m_currentTable[positionNewTable] == EMPTY) {

                        m_currentTable[positionNewTable] = m_oldTable[m_location];
                        m_oldTable[m_location] = DELETED;
                        ++m_currentSize;
                    }
                    else {

                        int count = 0;
                        bool isAdded = false;
                        int newIndexHash = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;

                        while (!isAdded) {

                            positionNewTable = (newIndexHash + (count * count)) % m_currentCap;
                            if (m_currentTable[positionNewTable] == EMPTY) {

                                m_currentTable[positionNewTable] = m_oldTable[m_location];
                                m_oldTable[m_location] = DELETED;
                                ++m_currentSize;
                                isAdded = true;
                            }
                            ++count;
                        }
                    }
                }
            }

            ++m_quarterChecker;
            break;

        case 2:
            //int positionNewTable = 0;

            for (; m_location < m_quarterSize * m_quarterChecker; m_location++) {

                if (m_oldTable[m_location].getSequence() != "" && m_oldTable[m_location].getSequence() != DELETEDKEY) {

                    positionNewTable = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;
                    if (m_currentTable[positionNewTable] == EMPTY) {

                        m_currentTable[positionNewTable] = m_oldTable[m_location];
                        m_oldTable[m_location] = DELETED;
                        ++m_currentSize;
                    }
                    else {

                        int count = 0;
                        bool isAdded = false;
                        int newIndexHash = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;

                        while (!isAdded) {

                            positionNewTable = (newIndexHash + (count * count)) % m_currentCap;
                            if (m_currentTable[positionNewTable] == EMPTY) {

                                m_currentTable[positionNewTable] = m_oldTable[m_location];
                                m_oldTable[m_location] = DELETED;
                                ++m_currentSize;
                                isAdded = true;
                            }
                            ++count;
                        }
                    }
                }
            }

            ++m_quarterChecker;
            break;

        case 3:
            //int positionNewTable = 0;

            for (; m_location < m_quarterSize + (m_quarterChecker * 3); m_location++) {

                if (m_oldTable[m_location].getSequence() != "" && m_oldTable[m_location].getSequence() != DELETEDKEY) {

                    positionNewTable = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;
                    if (m_currentTable[positionNewTable] == EMPTY) {

                        m_currentTable[positionNewTable] = m_oldTable[m_location];
                        m_oldTable[m_location] = DELETED;
                        ++m_currentSize;
                    }
                    else {

                        int count = 0;
                        bool isAdded = false;
                        int newIndexHash = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;

                        while (!isAdded) {

                            positionNewTable = (newIndexHash + (count * count)) % m_currentCap;
                            if (m_currentTable[positionNewTable] == EMPTY) {

                                m_currentTable[positionNewTable] = m_oldTable[m_location];
                                m_oldTable[m_location] = DELETED;
                                ++m_currentSize;
                                isAdded = true;
                            }
                            ++count;
                        }
                    }
                }
            }

            ++m_quarterChecker;
            break;

        case 4:
            //int positionNewTable = 0;

            for (; m_location < m_oldCap; m_location++) {

                if (m_oldTable[m_location].getSequence() != "" && m_oldTable[m_location].getSequence() != DELETEDKEY) {

                    positionNewTable = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;
                    if (m_currentTable[positionNewTable] == EMPTY) {

                        m_currentTable[positionNewTable] = m_oldTable[m_location];
                        m_oldTable[m_location] = DELETED;
                        ++m_currentSize;
                    }
                    else {

                        int count = 0;
                        bool isAdded = false;
                        int newIndexHash = m_hash(m_oldTable[m_location].getSequence()) % m_currentCap;

                        while (!isAdded) {

                            positionNewTable = (newIndexHash + (count * count)) % m_currentCap;
                            if (m_currentTable[positionNewTable] == EMPTY) {

                                m_currentTable[positionNewTable] = m_oldTable[m_location];
                                m_oldTable[m_location] = DELETED;
                                ++m_currentSize;
                                isAdded = true;
                            }
                            ++count;
                        }
                    }
                }
            }

            delete[] m_oldTable;
            m_quarterChecker = 0;
            break;

        default:
            break;
        }

    }

    //m_currentTable = new DNA(dna);
    if (isInserted) {
        return true;
    }
    return false;
}

bool DnaDb::remove(DNA dna) {

    return false;
}

DNA DnaDb::getDNA(string sequence, int location) {

    // quadratic probe to find the object
    if (m_isOldActivated) {

        int index = m_hash(sequence) % m_currentCap;
        if (m_currentTable[index].getSequence() == "" && m_oldTable[index].getSequence() ==  "") {

            return EMPTY;
        }

        if (m_currentTable[index].getSequence() == sequence && m_currentTable[index].getLocId() == location) {

            return m_currentTable[index];
        }

        else if (m_oldTable[index].getSequence() == sequence && m_oldTable[index].getLocId() == location) {

            return m_oldTable[index];
        }

        else {

            int count = 0;
            bool isFinished = false;
            int newIndexHash = m_hash(sequence) % m_currentCap;

            while (!isFinished) {

                index = (newIndexHash + (count * count)) % m_currentCap;
                if (m_currentTable[index].getSequence() == "" && m_oldTable[index].getSequence() == "") {

                    isFinished = true;
                    return EMPTY;
                }
                else {

                    if (m_currentTable[index].getSequence() == sequence && m_currentTable[index].getLocId() == location) {

                        isFinished = true;
                        return m_currentTable[index];
                    }

                    if (m_oldTable[index].getSequence() == sequence && m_oldTable[index].getLocId() == location) {

                        isFinished = true;
                        return m_oldTable[index];
                    }

                    ++count;
                }
            }
        }
    }
    else {

        int index = m_hash(sequence) % m_currentCap;
        if (m_currentTable[index].getSequence() == "") {

            return EMPTY;
        }

        if (m_currentTable[index].getSequence() == sequence && m_currentTable[index].getLocId() == location) {

            return m_currentTable[index];
        }

        else {

            int count = 0;
            bool isFinished = false;
            int newIndexHash = m_hash(sequence) % m_currentCap;

            while (!isFinished) {

                index = (newIndexHash + (count * count)) % m_currentCap;
                if (m_currentTable[index].getSequence() == "") {

                    isFinished = true;
                    return EMPTY;
                }
                else {

                    if (m_currentTable[index].getSequence() == sequence && m_currentTable[index].getLocId() == location) {

                        isFinished = true;
                        return m_currentTable[index];
                    }

                    ++count;
                }
            }
        }
    }
}

float DnaDb::lambda() const {

    return (m_currentSize / m_currentCap);
}

float DnaDb::deletedRatio() const {

    return (m_currNumDeleted / m_currentCap);
}

void DnaDb::dump() const {
    cout << "Dump for current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number) {
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current) {
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME - 1;
    for (int i = current; i < MAXPRIME; i++) {
        for (int j = 2; j * j <= i; j++) {
            if (i % j == 0)
                break;
            else if (j + 1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

DNA::DNA(string sequence, int location) {
    if ((location >= MINLOCID && location <= MAXLOCID) ||
        (location == 0 && sequence == "DELETED")) {
        // this is a normal or a DELETED object
        m_sequence = sequence;
        m_location = location;
    }
    else {
        // this is the empty object
        m_sequence = "";
        m_location = 0;
    }
}

string DNA::getSequence() const {
    return m_sequence;
}

int DNA::getLocId() const {
    return m_location;
}

// Overloaded assignment operator
const DNA& DNA::operator=(const DNA& rhs) {
    if (this != &rhs) {
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

// Overloaded insertion operator.  Prints DNA's sequence (key),
// and the location ID. This is a friend function in DNA class.
ostream& operator<<(ostream& sout, const DNA& dna) {
    if (!dna.m_sequence.empty())
        sout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
    else
        sout << "";
    return sout;
}

// Overloaded equality operator. This is a friend function in DNA class.
// To test inequality we may negate the results of this operator.
bool operator==(const DNA& lhs, const DNA& rhs) {
    return ((lhs.m_sequence == rhs.m_sequence) && (lhs.m_location == rhs.m_location));
}

bool DnaDb::isValidLocation(DNA dna) {

    if (dna.getLocId() >= MINLOCID && dna.getLocId() <= MAXLOCID) {

        return true;
    }
    else {
        return false;
    }
}