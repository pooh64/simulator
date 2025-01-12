#include "mapgenerator.h"
#include "locus.h"
#include "relaxationtor.h"
#include <QDebug>
#include <QRandomGenerator>
#include <QTimer>

MapGenerator::MapGenerator(uint32_t xSize, uint32_t ySize, QObject *parent)
    : RelaxationTor(xSize, ySize, parent) {
  m_hTimer = new QTimer(this);
  connect(m_hTimer, SIGNAL(timeout()), this, SLOT(heightRelaxation()));
  m_rTimer = new QTimer(this);
  connect(m_rTimer, SIGNAL(timeout()), this, SLOT(riverFlow()));
}

void MapGenerator::startRelaxation(int32_t repeat) {
  m_timer->stop();
  generateLocuses();
  m_repeat = repeat;
  m_timer->start(100);
}

void MapGenerator::stopRelaxation() {
  m_timer->stop();

  for (Locus &locus : m_locuses) {
    locus.clear();
  }
  gatherPoints();

  for (Locus &locus : m_locuses) {
    for (uint32_t &point : locus.m_points) {
      m_locusMap[point] = &locus;
    }
  }
  for (Locus &locus : m_locuses) {
    findNeighbors(locus);
  }
  drawSpace();
  emit sendSpace(m_space, m_xSize, m_ySize);
  emit finishRelax();
}

void MapGenerator::findNeighbors(Locus &locus) {
  for (uint32_t &point : locus.m_points) {
    uint32_t x = point % m_xSize;
    uint32_t y = point / m_xSize;
    uint32_t x_l = (x == 0) ? m_xSize - 1 : x - 1;
    uint32_t x_r = (x == m_xSize - 1) ? 0 : x + 1;
    uint32_t y_u = (y == 0) ? m_ySize - 1 : y - 1;
    uint32_t y_d = (y == m_ySize - 1) ? 0 : y + 1;
    locus.addNeighbor(m_locusMap[y_u * m_xSize + x_l]);
    locus.addNeighbor(m_locusMap[y_u * m_xSize + x]);
    locus.addNeighbor(m_locusMap[y_u * m_xSize + x_r]);
    locus.addNeighbor(m_locusMap[y * m_xSize + x_l]);
    locus.addNeighbor(m_locusMap[y * m_xSize + x_r]);
    locus.addNeighbor(m_locusMap[y_d * m_xSize + x_l]);
    locus.addNeighbor(m_locusMap[y_d * m_xSize + x]);
    locus.addNeighbor(m_locusMap[y_d * m_xSize + x_r]);
  }
}

void MapGenerator::generateHeight(int32_t repeat) {
  // Mountains
  for (Locus &locus : m_locuses) {
    if ((0.1 * m_xSize < locus.m_x) && (locus.m_x < 0.9 * m_xSize) &&
        (0.1 * m_ySize < locus.m_y) && (locus.m_y < 0.9 * m_ySize) &&
        QRandomGenerator::global()->bounded(m_locuses.size() / 30.0) == 0) {
      locus.m_z = 0x80 + (QRandomGenerator::global()->generate() & 0x7f);
      locus.m_fixZ = true;
    } else {
      locus.m_z = QRandomGenerator::global()->generate() & 0xff;
      locus.m_fixZ = false;
    }
  }
  // Water barrier
  for (uint32_t y = 0; y < m_ySize; y++) {
    m_locusMap[y * m_xSize]->m_z = 0;
    m_locusMap[y * m_xSize]->m_fixZ = true;
    m_locusMap[y * m_xSize + m_xSize - 1]->m_z = 0;
    m_locusMap[y * m_xSize + m_xSize - 1]->m_fixZ = true;
  }
  for (uint32_t x = 0; x < m_xSize; x++) {
    m_locusMap[x]->m_z = 0;
    m_locusMap[x]->m_fixZ = true;
    m_locusMap[(m_ySize - 1) * m_xSize + x]->m_z = 0;
    m_locusMap[(m_ySize - 1) * m_xSize + x]->m_fixZ = true;
  }
  drawSpace();
  emit sendSpace(m_space, m_xSize, m_ySize);
  m_repeat = repeat;
  m_hTimer->start(100);
}

void MapGenerator::heightRelaxation() {
  for (Locus &locus : m_locuses) {
    locus.averageZ();
  }
  for (Locus &locus : m_locuses) {
    locus.calcHeightColor();
  }
  if (m_repeat != -1) {
    m_repeat--;
    if (m_repeat == 0) {
      m_hTimer->stop();
      addMapRand();
    }
  }
  drawSpace();
  emit sendSpace(m_space, m_xSize, m_ySize);
}

void MapGenerator::addMapRand() {
  for (Locus &locus : m_locuses) {
    if (!locus.m_fixZ) {
      locus.m_z += QRandomGenerator::global()->generate() & 0xf;
      locus.calcHeightColor();
    }
  }
}

void MapGenerator::select(uint32_t x, uint32_t y) {
  if (m_locusMap.empty()) {
    return;
  }
  Locus *locus = m_locusMap[y * m_xSize + x];
  locus->m_color = 0xffff0000;
  locus->drawSpace(m_space);
  for (Locus *neighbor : locus->m_neighbors) {
    neighbor->m_color = 0xffffa000;
    neighbor->drawSpace(m_space);
  }
  emit sendSpace(m_space, m_xSize, m_ySize);
}

void MapGenerator::riverGeneration() {
  // River source
  for (Locus &locus : m_locuses) {
    if (!locus.m_fixZ && locus.m_z > 0x40 &&
        QRandomGenerator::global()->bounded(m_locuses.size() / 50.0) == 0) {
      locus.m_river = true;
    }
    locus.calcHeightColor();
  }
  drawSpace();
  emit sendSpace(m_space, m_xSize, m_ySize);
  m_rTimer->start(1000);
}

void MapGenerator::riverFlow() {
  std::vector<Locus *> newRiver;
  for (Locus &locus : m_locuses) {
    if (locus.m_river) {
      Locus *flow = locus.minNeighbor();
      if (!flow->m_river) {
        newRiver.push_back(flow);
      } else {
        if (flow->m_z > locus.m_z) {
          flow = locus.minNonRiver();
          if (flow != nullptr) {
            newRiver.push_back(flow);
          }
        }
      }
    }
  }
  if (newRiver.empty()) {
    m_rTimer->stop();
  }
  for (Locus *locus : newRiver) {
    locus->m_river = true;
    locus->calcHeightColor();
  }
  drawSpace();
  emit sendSpace(m_space, m_xSize, m_ySize);
}
