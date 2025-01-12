#ifndef CPU_H
#define CPU_H

#include "instr.h"
#include <QObject>
#include <QThread>

constexpr uint32_t REG_SIZE = 16;
constexpr uint32_t MEM_SIZE = 65536;
constexpr uint32_t DIS_WIDTH = 256;
constexpr uint32_t DIS_HEIGHT = 128;
constexpr uint32_t DIS_SCALE = 4;

class CPU;

class CPUExecThreads : public QThread {
public:
  explicit CPUExecThreads(CPU *cpu_init) { cpu = cpu_init; };
  CPU *cpu;
  void run();
};

class CPU : public QObject {
  Q_OBJECT
public:
  explicit CPU(QObject *parent = nullptr);
  ~CPU();
  std::string processLabels(std::string input, std::stringstream &labels);
  uint32_t &readMem(uint32_t addr);
public slots:
  void readInstrs(QString input);
  void run();
  void pause();
  void stop();
  void step();
  void dumpStatus();
  void dumpMem();
  void updateDisplay();
  void showMsg(QString msg);
signals:
  void statusUpd(QString status);
  void memUpd(QString mem);
  void sendMsg(QString msg);
  void displayUpd(uint32_t *mem, uint32_t width, uint32_t height,
                  uint32_t scale);

public:
  uint32_t m_regFile[REG_SIZE] = {};
  uint32_t m_mem[MEM_SIZE] = {};
  uint32_t m_PC;
  uint32_t m_nextPC;
  uint32_t m_run;
  CPUExecThreads *execThread;
};

#endif // CPU_H
