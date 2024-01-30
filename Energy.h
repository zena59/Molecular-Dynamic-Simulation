#pragma once
//NOTE : IN ALL CALCULATIONS,SIGMA IS TAKEN AS EQUAL TO UNITY


//force calulator
float force(float r)
{
	//eps is epsellon, well depth
	const float eps = 1;

	float F = 0;
	if (r >= 3)
	{
		F = 0;
	}
	else

	{
		float f;
		f = 4 * eps * ((12 / (pow(r, 13))) - (6 / (pow(r, 7))));

		F = f + 0.01094383;
	}
	return F;
}

//potential
float psi(float r)
{
	float potential = 0;
	if (r >= 3)
	{
		potential = 0;
	}

	else
	{
		const float eps = 1;
		float v;
		float calc = 1 / r;
		v = 4 * eps * (pow(calc, 12) - pow(calc, 6));
		potential = v + (0.01094383 * r) - 0.027352048;
	}
	return potential;
}