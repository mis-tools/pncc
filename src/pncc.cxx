#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm> 
#include <string> 
#include <cstring> 

using namespace std;

struct NormalizedCorrelation
{
    char m_acAtlasFilename[128];
    int m_iIndex;
    float m_fResult;
};

bool compare_lt( const NormalizedCorrelation &kA, const NormalizedCorrelation &kB)
{
    return kA.m_fResult < kB.m_fResult;
}
bool compare_gt( const NormalizedCorrelation &kA, const NormalizedCorrelation &kB)
{
    return kA.m_fResult > kB.m_fResult;
}

int main(const int argc, const char * argv[])
{
    if ( argc < 6 )
    {
        cout << "Usage: " << argv[0] << " InputImage.raw NumberOfImagesWanted AtlasBaseDirectory AtlasImageFilename AtlasID1 [AtlasID2 ... ]" << std::endl;
		cout << "       All images must have same resolution and reside in the same physical space." << std::endl;
		cout << "       Prints out ID's of the {NumberOfImagesWanted} atlas images most similar to the InputImage." << std::endl;
        return EXIT_FAILURE;
    }
    string kInputImage = argv[1];
    const int iNumberOfImagesWanted = atoi(argv[2]);
    string kAtlasBaseDirectory = argv[3];
    string kAtlasImgFilename = argv[4];  // image in mhd/raw format

    // Read image
    FILE * pkFile;
    pkFile = fopen(kInputImage.c_str(), "r");
    fseek(pkFile, 0, SEEK_END);
    const long lFileByteSize = ftell(pkFile);
    rewind(pkFile);
    float * pfImage = (float*)malloc(lFileByteSize);
    size_t iStat = fread(pfImage, 1, lFileByteSize, pkFile);
    fclose(pkFile);

    const long numElem = lFileByteSize / sizeof(float);
    const float fInvNumElem = 1.0f / float(numElem);

    vector<NormalizedCorrelation> kAtlasCorrelations;
    
    kAtlasCorrelations.resize(argc - 5);

    // Read atlas
    #pragma omp parallel for private(pkFile)
    for ( int n = 5; n < argc; n++ )
    {
        bool bUShortType = false; // Only float or ushort supported
        string kAtlasId = argv[n];
        string kAtlasHeaderFilename = kAtlasBaseDirectory + string("/") + kAtlasId + string("/") + kAtlasImgFilename;
        char acBuf[64];
        char acType[16];
        char acOrientation[3];
        int aiDimSize[3];
        FILE * pkHeaderFile = fopen(kAtlasHeaderFilename.c_str(), "r");
        char kAtlasImgFilenameRAW[256];
        while ( ! feof(pkHeaderFile) )
        {
            fgets(acBuf, 64, pkHeaderFile);
            if ( sscanf(acBuf, "DimSize = %d %d %d", &aiDimSize[0], &aiDimSize[1], &aiDimSize[2]) )
            {
                //printf("DimSize: %d %d %d\n", aiDimSize[0], aiDimSize[1], aiDimSize[2]);
            }
            if ( sscanf(acBuf, "ElementDataFile = %s", kAtlasImgFilenameRAW) )
            {
                //printf("ElementDataFile: %s\n", kAtlasImgFilenameRAW);
            }
            if ( sscanf(acBuf, "ElementType = %s", acType) )
            {
                if ( strcmp("MET_USHORT", acType) == 0 )
                {
                    bUShortType = true;
                }
            }
            if ( sscanf(acBuf, "AnatomicalOrientation = %s", acOrientation) )
            {
                if ( strcmp("RAI", acOrientation) != 0 ) {
                    cout << "wrong orientation of: " << kAtlasHeaderFilename;
                    cout << " = " << acOrientation;
                    cout << ", should have been = RAI" << endl;
                    exit(-1);
                }
            }
        }
        fclose(pkHeaderFile);

        const string kAtlasFilename = kAtlasBaseDirectory + string("/") + kAtlasId + string("/") + kAtlasImgFilenameRAW;
        pkFile = fopen(kAtlasFilename.c_str(), "r");
        float * pfAtlasImage;
        unsigned short * psAtlasImage;
        long alNumElem = aiDimSize[0] * aiDimSize[1] * aiDimSize[2];

        if (alNumElem != numElem) {
            cout << "wrong number of elements: " << alNumElem;
            cout << " != " << numElem;
            cout << endl;
            exit(-2);
        }

        if ( bUShortType )
        {
            psAtlasImage = (unsigned short *)malloc( alNumElem * sizeof(unsigned short));
            iStat = fread(psAtlasImage, alNumElem, sizeof(unsigned short), pkFile);
        }
        else
        {
            pfAtlasImage = (float *)malloc( alNumElem * sizeof(float));
            iStat = fread(pfAtlasImage, alNumElem, sizeof(float), pkFile);
        }
        fclose(pkFile);

        // calculate normalized cross correlation
        float fABProdSum = 0.0f;
        float fASum = 0.0f;
        float fBSum = 0.0f;
        float fASquaredSum = 0.0f;
        float fBSquaredSum = 0.0f;
        for ( int i = 0; i < alNumElem; i++ )
        {
            float fA = pfImage[i];
            float fB;
            if ( bUShortType )
            {
                  fB = static_cast<float>(psAtlasImage[i]);
            }
            else
            {
                fB = pfAtlasImage[i];
            }

            fABProdSum += fA * fB;
            fASum += fA;
            fBSum += fB;
            fASquaredSum += fA * fA;
            fBSquaredSum += fB * fB;
        }
        const float mean_s = fBSum * fInvNumElem;
        const float mean_t = fASum * fInvNumElem;
        const float var_s = fBSquaredSum * fInvNumElem - mean_s * mean_s;
        const float var_t = fASquaredSum * fInvNumElem - mean_t * mean_t;
        const float covariance = fABProdSum * fInvNumElem - mean_s*mean_t;

        float fResult;

        if ((var_s < 0.00001f) || (var_t < 0.00001f) )
        {
            fResult = -1.0f;
        }
        else
        {
            fResult = covariance / sqrt( var_s*var_t );
        }

        // Save result in struct
        strcpy(kAtlasCorrelations[n-5].m_acAtlasFilename, argv[n]);
        kAtlasCorrelations[n-5].m_iIndex = n;
        kAtlasCorrelations[n-5].m_fResult = fResult;
        if ( bUShortType )
        {
            free(psAtlasImage);
        }
        else
        {
            free(pfAtlasImage);
        }
    }
    #pragma omp barrier

    // Use the greater than compare function in order to reverse sort and
    // avoid reverse iteration.
    sort(kAtlasCorrelations.begin(), kAtlasCorrelations.end(), compare_gt);

    vector<NormalizedCorrelation>::const_iterator kIter;
    if (iNumberOfImagesWanted > 0) {
    for ( kIter = kAtlasCorrelations.begin(); kIter - kAtlasCorrelations.begin() < iNumberOfImagesWanted; kIter++ )
    {
        cout << (*kIter).m_acAtlasFilename << " " << flush;
    }
    cout << endl;
    } else {  // if iNumberOfImagesWanted is 0 the output values for for all images
        for (kIter = kAtlasCorrelations.begin(); kIter != kAtlasCorrelations.end(); kIter++) {
            cout << (*kIter).m_fResult << " " <<(*kIter).m_acAtlasFilename << endl;
        }
    }

    free(pfImage);
    return EXIT_SUCCESS;
}
