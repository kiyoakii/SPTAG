// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#include "inc/Core/VectorSet.h"
#include "inttypes.h"
#include <fstream>
#include <memory>
#include "inc/Helper/VectorSetReader.h"
#include "inc/Core/Common/CommonUtils.h"

using namespace SPTAG;

#pragma warning(disable:4996)  // 'fopen': This function or variable may be unsafe. Consider using fopen_s instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. See online help for details.

VectorSet::VectorSet()
{
}


VectorSet::~VectorSet()
{
}


BasicVectorSet::BasicVectorSet(const ByteArray& p_bytesArray,
                               VectorValueType p_valueType,
                               DimensionType p_dimension,
                               SizeType p_vectorCount)
    : m_data(p_bytesArray),
      m_valueType(p_valueType),
      m_dimension(p_dimension),
      m_vectorCount(p_vectorCount),
      m_perVectorDataSize(static_cast<SizeType>(p_dimension * GetValueTypeSize(p_valueType)))
{
}

void BasicVectorSet::readXvec(const char* p_filePath, VectorValueType p_valueType,
    DimensionType p_dimension, SizeType p_vectorCount) 
{
    size_t vectorDataSize = GetValueTypeSize(p_valueType) * p_dimension;
    size_t totalRecordVectorBytes = vectorDataSize * p_vectorCount;
    ByteArray l_data = std::move(ByteArray::Alloc(totalRecordVectorBytes));
    char* vecBuf = reinterpret_cast<char*>(l_data.Data());

    std::ifstream in(p_filePath, std::ifstream::binary);
    if (!in.is_open()) {
        fprintf(stderr, "Error: Failed to read input file: %s \n", p_filePath);
        exit(-1);
    }

    DimensionType dim = p_dimension;
    for (size_t i = 0; i < p_vectorCount; i++) {
        in.read((char*)&dim, 4);
        if (dim != p_dimension) {
            fprintf(stderr, "Error: Xvec file %s has No.%zd vector whose dims are not as many as expected. Expected: %d, Fact: %d\n", p_filePath, i, p_dimension, dim);
            exit(-1);
        }
        in.read(vecBuf + i * vectorDataSize, vectorDataSize);
    }

    in.close();

    m_data = std::move(l_data);
    m_valueType = p_valueType;
    m_dimension = p_dimension;
    m_vectorCount = p_vectorCount;
    m_perVectorDataSize = vectorDataSize;
}

void BasicVectorSet::readDefault(const char* p_filePath, VectorValueType p_valueType) 
{
    std::ifstream in(p_filePath, std::ifstream::binary);
    if (!in.is_open()) {
        fprintf(stderr, "Error: Failed to read input file: %s \n", p_filePath);
        exit(-1);
    }

    SizeType row;
    DimensionType col;
    in.read((char*)&row, 4);
    in.read((char*)&col, 4);

    size_t vectorDataSize = GetValueTypeSize(p_valueType) * col;
    size_t totalRecordVectorBytes = vectorDataSize * row;
    ByteArray l_data = std::move(ByteArray::Alloc(totalRecordVectorBytes));
    char* vecBuf = reinterpret_cast<char*>(l_data.Data());
    in.read(vecBuf, totalRecordVectorBytes);

    in.close();

    m_data = std::move(l_data);
    m_valueType = p_valueType;
    m_dimension = col;
    m_vectorCount = row;
    m_perVectorDataSize = vectorDataSize;
}

void BasicVectorSet::readTxt(const char* p_filePath, VectorValueType p_valueType, DimensionType p_dimension, std::string p_delimiter) {
    std::shared_ptr<SPTAG::Helper::ReaderOptions> options = std::make_shared<SPTAG::Helper::ReaderOptions>(p_valueType, p_dimension, p_delimiter, 1);
    auto vectorReader = SPTAG::Helper::VectorSetReader::CreateInstance(options);
    if (ErrorCode::Success != vectorReader->LoadFile(p_filePath))
    {
        fprintf(stderr, "Failed to read input file.\n");
        exit(1);
    }

	std::shared_ptr<VectorSet> ptr = vectorReader->GetVectorSet();
    BasicVectorSet* vectors = (BasicVectorSet*)ptr.get();

	m_data = std::move(vectors->m_data);
    m_valueType = vectors->m_valueType;
    m_dimension = vectors->m_dimension;
    m_vectorCount = vectors->m_vectorCount;
    m_perVectorDataSize = vectors->m_perVectorDataSize;
}

// copied from src/IndexBuilder/main.cpp
BasicVectorSet::BasicVectorSet(const char* p_filePath, VectorValueType p_valueType,
    DimensionType p_dimension, SizeType p_vectorCount, VectorFileType p_fileType, std::string p_delimiter, DistCalcMethod p_distCalcMethod)
{
    if (p_fileType == VectorFileType::XVEC)
    {
        readXvec(p_filePath, p_valueType, p_dimension, p_vectorCount);
    }
    else if (p_fileType == VectorFileType::DEFAULT)
    {
        readDefault(p_filePath, p_valueType);
    }
    else if (p_fileType == VectorFileType::TXT) 
    {
        readTxt(p_filePath, p_valueType, p_dimension, p_delimiter);
    }
    else
    {
        fprintf(stderr, "VectorFileType Unsupported.\n");
        exit(-1);
    }

    if (p_distCalcMethod == DistCalcMethod::Cosine) {
#pragma omp parallel for
        for (int64_t i = 0; i < m_vectorCount; i++)
        {
            int64_t offset = i * m_perVectorDataSize;
            switch (p_valueType)
            {
#define DefineVectorValueType(Name, Type) \
case SPTAG::VectorValueType::Name: \
SPTAG::COMMON::Utils::Normalize<Type>(reinterpret_cast<Type *>(m_data.Data() + offset), p_dimension, SPTAG::COMMON::Utils::GetBase<Type>()); \
break; \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
            default:
                break;
            }
        }
    }
}

BasicVectorSet::~BasicVectorSet()
{
}


VectorValueType
BasicVectorSet::GetValueType() const
{
    return m_valueType;
}


void*
BasicVectorSet::GetVector(SizeType p_vectorID) const
{
    if (p_vectorID < 0 || p_vectorID >= m_vectorCount)
    {
        return nullptr;
    }

    return reinterpret_cast<void*>(m_data.Data() + ((size_t)p_vectorID) * m_perVectorDataSize);
}


void*
BasicVectorSet::GetData() const
{
    return reinterpret_cast<void*>(m_data.Data());
}

DimensionType
BasicVectorSet::Dimension() const
{
    return m_dimension;
}


SizeType
BasicVectorSet::Count() const
{
    return m_vectorCount;
}


bool
BasicVectorSet::Available() const
{
    return m_data.Data() != nullptr;
}


ErrorCode 
BasicVectorSet::Save(const std::string& p_vectorFile) const
{
    FILE * fp = fopen(p_vectorFile.c_str(), "wb");
    if (fp == NULL) return ErrorCode::FailedOpenFile;

    fwrite(&m_vectorCount, sizeof(SizeType), 1, fp);
    fwrite(&m_dimension, sizeof(DimensionType), 1, fp);

    fwrite((const void*)(m_data.Data()), m_data.Length(), 1, fp);
    fclose(fp);
    return ErrorCode::Success;
}

SizeType BasicVectorSet::PerVectorDataSize() const {
    return (SizeType)m_perVectorDataSize;
}