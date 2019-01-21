#include <iostream>

using namespace std;

class Cal_Pi
{
public:
	double TaylorPi(int k);	//kΪ����
	double BBPPi(int k);	//kΪ��ȷ����λ��������ȡk=15
};

class Cal_lnPi
{
	//�����Ϊ[a,b]��nΪ�ֶ���
public:
	double Compound_trapezium(double a, double b, int n);	
	double Compound_Simpson2(double a, double b, int n);	//ֱ���ù�ʽ��
	//�������㷨
	double Compound_Simpson(double a, double b, int n);
	double Cotes(double a, double b, int n);
	double Romberg(double a, double b, int n);
	double Newton(double x_unit, double pi, int nCount);
};

class Cal_Pix
{
	//�����Ϊ[a,b]��nΪ�ֶ���
public:
	//RK
	double Runge_Kutta4(double a, double b, int n);
	double CalPixRK(int n, double exp, double pi, double lnpi);
	//Taylor
	double CalPixTaylor(double exp, double pi, double lnpi);
};

//�������3 ţ�ٷ�
double Cal_sqrt3(int n)
{
	double result = 2;
	for (int i = 0; i < n; i++)
		result -= (result * result - 3) / (2 * result);
	return result;
}

//����e^x ̩��չ��
double Myexp(double x, int nCount)
{
	double result = 1;
	double Factorial = 1;
	double xinput = x;
	for (int i = 1; i <= nCount; i++)
	{
		Factorial *= i;
		result += 1.0 * xinput / Factorial;
		xinput *= x;
	}
	return result;
}

double Cal_Pi::TaylorPi(int k)
{
	double result = 0;
	double sqrt3 = Cal_sqrt3(6);	//��¼sqrt3/3
	double x = 1;
	for (int i = 0; i < k; i++)
	{
		if (i % 2 == 0)
		{
			result += 1.0 / x / (2 * i + 1);
		}
		else
		{
			result -= 1.0 / x / (2 * i + 1);
		}
		x *= 3;
	}
	return 2 * sqrt3 * result;
}

//k=15 kΪλ������
double Cal_Pi::BBPPi(int k)
{
	if (k < 0)
		k = 0;
	double result = 0;
	for (int i = 0; i < k; i++)
	{
		/*double temp = 4.0 / (8 * i + 1) - 2.0 / (8 * i + 4) -
			1.0 / (8 * i + 5) - 1.0 / (8 * i + 6);*/
		double temp = 1.0 * (120 * i * i + 151 * i + 47) / (512 * i*i*i*i +
			1024 * i*i*i + 712 * i*i + 194 * i + 15);
		for (int j = 0; j < i; j++)
		{
			temp /= 16.0;
		}
		result += temp;
	}
	return result;
}

//�������ι�ʽ	1/x����
double Cal_lnPi::Compound_trapezium(double a, double b, int n)
{
	double result = 0;
	for (int k = 0; k <= n; k++)
	{
		if (k == 0)
			result += 1.0 / a;
		else if (k == n)
			result += 1.0 / b;
		else
			result += 2.0 * 1.0 / (a + k * (b - a) / n);
	}
	return 0.5 * (b - a) / n * result;
}

//������������ʽ
double Cal_lnPi::Compound_Simpson2(double a, double b, int n)
{
	if (n <= 0)
		n = 1;
	double result = 0;
	for (int k = 0; k <= n; k++)
	{
		if (k == 0)
			result += 1.0 / a + 4.0 / (a + k * (b - a) / n + 0.5 * (b - a) / n);
		else if (k == n)
			result += 1.0 / b;
		else
			result += 2.0 / (a + k * (b - a) / n) + 4.0 / (a + k * (b - a) / n + 0.5 * (b - a) / n);
	}
	return 1.0 / 6 * (b - a) / n * result;
}

//������������ʽ���������㷨��
double Cal_lnPi::Compound_Simpson(double a,double b,int n)
{
	return 4.0 / 3 * Compound_trapezium(a, b, 2 * n) - 1.0 / 3 * Compound_trapezium(a, b, n);
}

//����˹��ʽ
double Cal_lnPi::Cotes(double a, double b, int n)
{
	return 16.0 / 15 * Compound_Simpson(a, b, 2 * n) - 1.0 / 15 * Compound_Simpson(a, b, n);
}

//������ʽ
double Cal_lnPi::Romberg(double a, double b, int n)
{
	return 64.0 / 63 * Cotes(a, b, n) - 1.0 / 63 * Cotes(a, b, n);
}

//ţ�ٷ� e^x = pi
double Cal_lnPi::Newton(double x_unit, double pi, int nCount)
{
	double result = x_unit;
	for (int i = 0; i < nCount; i++)
		result += pi / Myexp(result,21) - 1;
	return result;
}


//�Ľ����������ʽ
//y(0) = 1,y_{n+1} = (1 + h + h^2/2 + h^3/6 + h^4/24)y_n
double Cal_Pix::Runge_Kutta4(double a, double b, int n)
{
	double result = 1;	//y(0) = 1
	double h = 1.0 * (b - a) / n;
	for (int i = 0; i < n; i++)
		result *= (1 + h + h*h / 2 + h*h*h / 6 + h*h*h*h / 24);
	return result;
}

double Cal_Pix::CalPixRK(int n, double exp, double pi, double lnpi)
{
	int expi = (int)exp;
	double expj = exp - expi;
	double rslttmp = Runge_Kutta4(0, expj*lnpi, n);
	double result = 1;
	for (int i = 0; i < expi; i++)
		result *= pi;
	return result * rslttmp;
}

double Cal_Pix::CalPixTaylor(double exp, double pi, double lnpi)
{
	int expi = (int)exp;
	double expj = exp - expi;
	double rslttmp = Myexp(expj * lnpi, 20);
	double result = 1;
	for (int i = 0; i < expi; i++)
		result *= pi;
	return result * rslttmp;
}

int main()
{
	Cal_Pi SoluPi;
	Cal_lnPi SoluLnPi;
	Cal_Pix SoluPix;
	printf("��Pi��\n");
	double pi = SoluPi.TaylorPi(30);
	printf("Taylor_pi��%.16f\n", pi);
	pi = SoluPi.BBPPi(12);
	printf("BBP_pi��%.16f\n", pi);

	printf("��Ln(pi)\n");
	double Lnx = SoluLnPi.Compound_Simpson2(1, pi, 8192);
	printf("Compound_Simpson��%.16f\n", Lnx);
	Lnx = SoluLnPi.Cotes(1, pi, 256);
	printf("Cotes��%.16f\n", Lnx);
	Lnx = SoluLnPi.Romberg(1, pi, 128);
	printf("Romberg��%.16f\n", Lnx);
	Lnx = SoluLnPi.Newton(1.5, pi, 5);
	printf("Newton��%.16f\n", Lnx);

	double x = 0;
	while (true)
	{
		printf("��pi^x\n");
		printf("������x��ֵ[1,10]������-1�˳���\n");
		cin >> x;
		if (x == -1)
			break;
		x = x > 10 ? 10 : x;
		x = x < 0 ? 0 : x;
		double Pix = SoluPix.CalPixRK(1024, x, pi, Lnx);
		printf("RK4��%.14f\n", Pix);
		Pix = SoluPix.CalPixTaylor(x, pi, Lnx);
		printf("Taylor��%.14f\n\n", Pix);
		//printf("%.14f\n", SoluPix.Runge_Kutta4(0, x*Lnx, 8192));
	}
	return 0;

}



