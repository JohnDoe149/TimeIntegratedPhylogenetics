#include "Alignment.hpp"
#include "Msg.hpp"
#include "ncl/nxscharactersblock.cpp"

Alignment::Alignment(std::string fn) {
    
    MultiFormatReader nexusReader;
    const char* fileName = fn.c_str();
    nexusReader.ReadFilepath(fileName, MultiFormatReader::NEXUS_FORMAT);

    size_t numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    if(numTaxaBlocks > 1)
        Msg::error("Too many taxa blocks (> 1)");

    for(size_t tBlock = 0; tBlock < numTaxaBlocks; tBlock++){
        NxsTaxaBlock* taxaBlock = nexusReader.GetTaxaBlock(tBlock);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        const unsigned numCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        const unsigned numUnalignedCharBlocks = nexusReader.GetNumUnalignedBlocks(taxaBlock);
        
        if(numUnalignedCharBlocks > 0)
            Msg::error("No unaligned data allowed!");
        if(numCharBlocks > 1)
            Msg::error("Too many char blocks (> 1)");

        for(size_t cBlock = 0; cBlock < numCharBlocks; cBlock++){
            NxsCharactersBlock* charBlock = nexusReader.GetCharactersBlock(taxaBlock, cBlock);
            std::string charBlockTitle = charBlock->GetTitle();
            stateSpace = 4;
            for(int i = 0; i < stateSpace; i++)
                frequencies.push_back(0);
            readNucleotideData(charBlock);
            std::cout << "Using nucleotide data type." << std::endl;
        }
        
    }
}

Alignment::~Alignment(){
    delete [] matrix[0];
    delete [] matrix;
}

Matrix<double> Alignment::getPairwiseIdentities(){
    Matrix<double> pairwiseID(numTaxa, 0.0);

    for(int i = 0; i < numTaxa; i++){
        for(int j = i; j < numTaxa; j++){
            if(i == j)
                pairwiseID(i, j) = 1.0;
            else{
                double pID = 0.0;
                for(int k = 0; k < numChar; k++)
                    pID += (matrix[i][k] & matrix[j][k]) > 0; // This implementation isn't ambiguous-character-friendly
                
                pID = pID/numChar;
                pairwiseID(i, j) = pID;
                pairwiseID(j, i) = pID;
            }
        }
    }

    return pairwiseID;
}

void Alignment::readNucleotideData(NxsCharactersBlock* charBlock){
    numTaxa = charBlock->GetNumActiveTaxa();
    numChar = charBlock->GetNumActiveChar();

    matrix = new unsigned long long int*[numTaxa];
    matrix[0] = new unsigned long long int[numTaxa*numChar];
    for(int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1] + numChar;
    for(int i = 0; i < numTaxa; i++)
        for(int j = 0; j < numChar; j++)
            matrix[i][j] = 0;

    for(int i = 0; i < numTaxa; i++){
        taxaNames.push_back(charBlock->GetTaxonLabel(i));
        for(int j = 0; j < numChar; j++){
            char state = charBlock->GetState(i, j);
            if(state == 'A' || state == 'R' || state == 'W' || state == 'M' || state == 'D' || state == 'H' || state == 'V'){
                matrix[i][j] += 1;
                frequencies[0] += 1;
            }
            if(state == 'C' || state == 'Y' || state == 'S' || state == 'M' || state == 'B' || state == 'H' || state == 'V'){
                matrix[i][j] += 2;
                frequencies[1] += 1;
            }
            if(state == 'G' || state == 'R' || state == 'S' || state == 'K' || state == 'B' || state == 'D' || state == 'V'){
                matrix[i][j] += 4;
                frequencies[2] += 1;
            }
            if(state == 'T' || state == 'Y' || state == 'W' || state == 'K' || state == 'B' || state == 'D' || state == 'H' || state == 'U'){
                matrix[i][j] += 8;
                frequencies[3] += 1;
            }
            if(state == 'N' || charBlock->IsGapState(i, j) || charBlock->IsMissingState(i, j))
                matrix[i][j] = 15;
        }
    }

    double total = 0.0;
    for(double v : frequencies){
        total += v;
    }

    for(double& v : frequencies){
        v = v/total;
    }
}