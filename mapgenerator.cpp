#include "mapgenerator.h"
#include "locus.h"
#include "relaxationtor.h"
#include <QDebug>

MapGenerator::MapGenerator(uint32_t xSize, uint32_t ySize, QObject *parent) :
    QObject(parent), m_xSize(xSize), m_ySize(ySize)
{
    m_space = new uint32_t[xSize*ySize];
}

MapGenerator::~MapGenerator() {
    delete m_space;
}

void MapGenerator::generateLocuses() {
    RelaxationTor relaxTor(m_xSize, m_ySize, this);
    relaxTor.generateLocuses();
    relaxTor.relax();
    relaxTor.relax();
    relaxTor.copyLocuses(&m_locuses);
    for (Locus& locus: m_locuses) {
        for (uint32_t& point: locus.m_points) {
            m_locusMap[point] = &locus;
        }
    }
    for (Locus& locus: m_locuses) {
        findNeighbors(locus);
    }
    drawSpace();
    emit sendSpace(m_space, m_xSize, m_ySize);
}

void MapGenerator::findNeighbors(Locus &locus) {
    for (uint32_t& point: locus.m_points) {
        uint32_t x = point % m_xSize;
        uint32_t y = point / m_xSize;
        uint32_t x_l = (x == 0) ? m_xSize - 1 : x - 1;
        uint32_t x_r = (x == m_xSize - 1) ? 0 : x + 1;
        uint32_t y_u = (y == 0) ? m_ySize - 1 : y - 1;
        uint32_t y_d = (y == m_ySize - 1) ? 0 : y + 1;
        locus.addNeighbor(m_locusMap[y_u*m_xSize + x_l]);
        locus.addNeighbor(m_locusMap[y_u*m_xSize + x_r]);
        locus.addNeighbor(m_locusMap[y_d*m_xSize + x_l]);
        locus.addNeighbor(m_locusMap[y_d*m_xSize + x_r]);
    }
}

void MapGenerator::generateHeight() {
    for(Locus& locus: m_locuses) {
        locus.m_z = qrand() & 0xff;
        locus.calcColor();
    }
    drawSpace();
    emit sendSpace(m_space, m_xSize, m_ySize);
}

void MapGenerator::select(uint32_t x, uint32_t y) {
    Locus* locus = m_locusMap[y*m_xSize + x];
    locus->m_color = 0xffff0000;
    locus->drawSpace(m_space);
    for (Locus* neighbor: locus->m_neighbors) {
        neighbor->m_color = 0xffffa000;
        neighbor->drawSpace(m_space);
    }
    emit sendSpace(m_space, m_xSize, m_ySize);
}

void MapGenerator::drawSpace() {
    for(Locus& locus: m_locuses) {
        locus.drawSpace(m_space);
    }
}
