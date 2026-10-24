# 学习增强 Treap 实验 —— 构建脚本
#   用法：make            （编译 la_treap.exe）
#         make clean      （清理产物）
#   说明：-static-libstdc++/-static-libgcc 静态链接运行时，
#         避免 PATH 中其他 MinGW 的 libstdc++-6.dll 造成 ABI 不匹配（见 runner.h）。

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -static-libstdc++ -static-libgcc
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
