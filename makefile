main: main.o ThomasSolver.o ImpliedVolatilitySurface.o DupireLocalVolatilitySurface.o Complex.o TrapezoidalRule.o HestonModel.o HestonPricer.o Chisquare.o HestonLocalVolatility.o RootSearch.o 
	g++ main.o ThomasSolver.o ImpliedVolatilitySurface.o DupireLocalVolatilitySurface.o Complex.o TrapezoidalRule.o HestonModel.o HestonPricer.o Chisquare.o HestonLocalVolatility.o RootSearch.o -o testing

main.o: main.cpp
	g++ -c main.cpp

ThomasSolver.o : ThomasSolver.cpp
	g++ -c ThomasSolver.cpp

ImpliedVolatilitySurface.o : ImpliedVolatilitySurface.cpp
	g++ -c ImpliedVolatilitySurface.cpp

DupireLocalVolatilitySurface.o : DupireLocalVolatilitySurface.cpp
	g++ -c DupireLocalVolatilitySurface.cpp

Complex.o : Complex.cpp
	g++ -c Complex.cpp

TrapezoidalRule.o : TrapezoidalRule.cpp
	g++ -c TrapezoidalRule.cpp

HestonModel.o : HestonModel.cpp
	g++ -c HestonModel.cpp

HestonPricer.o : HestonPricer.cpp
	g++ -c HestonPricer.cpp

Chisquare.o : Chisquare.cpp
	g++ -c Chisquare.cpp

HestonLocalVolatility.o : HestonLocalVolatility.cpp
	g++ -c HestonLocalVolatility.cpp

RootSearch.o : RootSearch.cpp	
	g++ -c RootSearch.cpp

clean:
	rm *.o testing


