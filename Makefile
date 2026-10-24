# 学习增强 Treap 实验 —— 构建脚本
#   用法：make            （编译 la_treap.exe）
#         make clean      （清理产物）
#   说明：-static 全静态链接 C++/pthread 运行时，exe 不再依赖 libstdc++-6.dll /
#         libwinpthread-1.dll，避免 PATH 中其他 MinGW 的同名 DLL 版本不一致导致段错误。

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -static
TARGET   := la_treap.exe
SRCS     := $(wildcard src/*.cpp)
OBJS     := $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean
