#include "OneWayANOVA.h"
#include "Factor.h"
#include "Response.h"

int main(int argc, char** argv)
{
	vector<double> dutpsa;

	dutpsa.push_back(0.003957);
	dutpsa.push_back(0.0038105);
	dutpsa.push_back(0.0036639);
	dutpsa.push_back(0.0033708);
	dutpsa.push_back(0.014216);

	dutpsa.push_back(0.013972);
	dutpsa.push_back(0.013776);
	dutpsa.push_back(0.014265);
	dutpsa.push_back(0.035369);
	dutpsa.push_back(0.035027);

	dutpsa.push_back(0.035173);
	dutpsa.push_back(0.034245);
	dutpsa.push_back(0.062531);
	dutpsa.push_back(0.06126);
	dutpsa.push_back(0.060088);

	dutpsa.push_back(0.053737);
	dutpsa.push_back(0.011724);
	dutpsa.push_back(0.012018);
	dutpsa.push_back(0.012262);
	dutpsa.push_back(0.011871);

	dutpsa.push_back(0.035711);
	dutpsa.push_back(0.039082);
	dutpsa.push_back(0.038935);
	dutpsa.push_back(0.041622);
	dutpsa.push_back(0.081094);

	dutpsa.push_back(0.078749);
	dutpsa.push_back(0.077675);
	dutpsa.push_back(0.06937);
	
	Response dutpsa_res(dutpsa);

	// Create BIAS Factor
	vector<int> bias;

	bias.push_back(0);
	bias.push_back(1);
	bias.push_back(2);
	bias.push_back(3);
	bias.push_back(3);

	bias.push_back(2);
	bias.push_back(1);
	bias.push_back(0);
	bias.push_back(0);
	bias.push_back(1);

	bias.push_back(2);
	bias.push_back(3);
	bias.push_back(3);
	bias.push_back(2);
	bias.push_back(1);

	bias.push_back(0);
	bias.push_back(0);
	bias.push_back(3);
	bias.push_back(2);
	bias.push_back(1);
	
	bias.push_back(0);
	bias.push_back(2);
	bias.push_back(1);
	bias.push_back(3);
	bias.push_back(3);

	bias.push_back(2);
	bias.push_back(1);
	bias.push_back(0);

	Factor bias_fac(bias, 4, Response(dutpsa));

	// Create DRG Factor
	vector<int> drg;

	drg.push_back(0);
	drg.push_back(0);
	drg.push_back(0);
	drg.push_back(0);
	drg.push_back(1);

	drg.push_back(1);
	drg.push_back(1);
	drg.push_back(1);
	drg.push_back(1);
	drg.push_back(1);

	drg.push_back(1);
	drg.push_back(1);
	drg.push_back(2);
	drg.push_back(2);
	drg.push_back(2);

	drg.push_back(2);
	drg.push_back(1);
	drg.push_back(1);
	drg.push_back(1);
	drg.push_back(1);

	drg.push_back(2);
	drg.push_back(2);
	drg.push_back(2);
	drg.push_back(2);
	drg.push_back(2);

	drg.push_back(2);
	drg.push_back(2);
	drg.push_back(2);

	Factor drg_fac(drg, 3, dutpsa_res);

	// Create Pulse Width factor
	vector<int> pw;

	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);

	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);

	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);
	pw.push_back(1);

	pw.push_back(1);
	pw.push_back(0);
	pw.push_back(0);
	pw.push_back(0);
	pw.push_back(0);

	pw.push_back(0);
	pw.push_back(0);
	pw.push_back(0);
	pw.push_back(0);
	pw.push_back(0);

	pw.push_back(0);
	pw.push_back(0);
	pw.push_back(0);
	
	Factor pw_fac(pw, 2, dutpsa_res);

	vector<Factor> facs;
	facs.push_back(bias_fac);
	facs.push_back(drg_fac);
	facs.push_back(pw_fac);

	OneWayANOVA mea(facs, dutpsa_res);

	std::cout << "For DUTPSA Data" << std::endl;
	for(int i = 0; i < 3 ; i++)
	{
		mea.printANOVATable(i);
		std::cout << std::endl;
	}
	return 0;
}
