/*
 *  RootFinding.c
 *  
 *
 *  Created by Garrett Weaver on 4/12/11.
 *
 */

#include <stdio.h>
#include <math.h>

#define	fpt double
#define currentFunct function_one
#define currentFunct_derivative function_one_derivative

//function_one is a polynomial that crosses though -3, 1 and 5
fpt function_one(fpt x)
{
	return ((x * x) - (4 * x) + 4.1) * (x + 1);
}
//function_one_derivative is the derivative function to 
fpt function_one_derivative(fpt x)
{
	return (3 * x * x) - (6 * x) + 0.1;
}
/*	this function is the bisection method 
 *	it returns the point where the function hits zero, 
 *	it needs 2 points on a function, one above and one below zero, denoted by x1 and x2
 *	it also needs a tolerence value that indicates how close we need to be to zero to be able to stop
 *	the second last value in the parameter set is a pointer to an int which indicates the number of iterations the function takes to give an answer
 *	the last value is to provide a stoping point for the loop incase we enter a case where we enter an infinite loop
 */
fpt bisectionMethod (fpt x1, fpt x2, fpt tolerence, int* iterations, int maxIterations)
{
	fpt xTemp1 = x1;
	fpt xTemp2 = x2;
	//xTemp3 is the middle point between x1 and x2
	fpt xTemp3 = (xTemp2 + xTemp1)/2;
	
	//set iterations equal to zero (one because we already made the calculation)
	*iterations = 1;
	
	//if the function crosses the 0 axis
	if((currentFunct(x1) * currentFunct(x2)) < 0)
	{
		//iterate until we are within a resonable tolerence
		while((fabs(currentFunct(xTemp3)) > tolerence) && (*iterations < maxIterations))
		{
			//if the root is between the mid point and the highest point
			if(currentFunct(xTemp3) * currentFunct(xTemp2) < 0)
			{
				//reset the lowest point to be the middle point
				xTemp1 = xTemp3;
			}
			//if the root is between the mid point and the lowest point
			else 
			{
				//reset the highest point to be the middle point
				xTemp2 = xTemp3;
			}
			
			//reset the middle point to be between the new highest and lowest points
			xTemp3 = (xTemp2 + xTemp1)/2;
			
			//add one to 'iterations'
			*iterations = *iterations + 1;
		}
		return xTemp3;
	}
	//if not, return a number out side of the range inputted to indicate that an error occured
	return (fabs(x1) + fabs(x2));
}

/*	this function is the regula-Falsi method
 *	it returns the point where the function hits zero, 
 *	it needs 2 points on a function, one above and one below zero, denoted by x1 and x2
 *	it also needs a tolerence value that indicates how close we need to be to zero to be able to stop
 *	the second last value in the parameter set is a pointer to an int which indicates the number of iterations the function takes to give an answer
 *	the last value is to provide a stoping point for the loop incase we enter a case where we enter an infinite loop
 */
fpt regulaFalsiMethod (fpt x1, fpt x2, fpt tolerence, int* iterations, int maxIterations)
{
	fpt xTemp1 = x1;
	fpt xTemp2 = x2;
	//xInter is the intersection point of the strait line where it crosses zero
	fpt xInter;
	
	//set iterations equal to zero
	*iterations = 0;
	
	//if the function crosses the 0 axis
	if((currentFunct(xTemp1) * currentFunct(xTemp2)) < 0)
	{
		do
		{
			//calculate the new x value
			xInter = xTemp1 + (currentFunct(xTemp1) * ((xTemp1 - xTemp2)/(currentFunct(xTemp2) - currentFunct(xTemp1))));
			
			//if the new x value's function return is below zero...
			if(currentFunct(xInter) < 0)
			{
				//if xTemp1 function's value is below zero
				if(currentFunct(xTemp1) < 0)
				{
					//replace xTemp1 with xInter
					xTemp1 = xInter;
				}
				//if xTemp2 function's value is below zero
				else 
				{
					//replace xTemp2 with xInter
					xTemp2 = xInter;
				}
			}
			//if the new x value's function return is above zero
			else 
			{
				//if xTemp1 function's value is above zero
				if(currentFunct(xTemp1) > 0)
				{
					//replace xTemp1 with xInter
					xTemp1 = xInter;
				}
				//if xTemp2 function's value is above zero
				else 
				{
					//replace xTemp2 with xInter
					xTemp2 = xInter;
				}
			}
			*iterations = *iterations + 1;
		}
		while((fabs(currentFunct(xInter)) > tolerence) && (*iterations < maxIterations));
		
		return xInter;
	}
	//if not, return a number outside of the range inputted to indicate that an error occured
	return (fabs(xTemp1) + fabs(xTemp2));
}

/*	this function is the secant method
 *	it returns the point where the function hits zero, 
 *	it needs 2 points on a function, one above and one below zero, denoted by x1 and x2
 *	it also needs a tolerence value that indicates how close we need to be to zero to be able to stop
 *	the second last value in the parameter set is a pointer to an int which indicates the number of iterations the function takes to give an answer
 *	the last value is to provide a stoping point for the loop incase we enter a case where we enter an infinite loop
 */
fpt secantMethod (fpt x1, fpt x2, fpt tolerence, int* iterations, int maxIterations)
{
	fpt xTemp1 = x1;
	fpt xTemp2 = x2;
	//xInter is the intersection point of the strait line where it crosses zero
	fpt xInter;
	
	//set iterations equal to zero
	*iterations = 0;
	
	do
	{
		//calculate the new x value
		xInter = xTemp1 + (currentFunct(xTemp1) * ((xTemp1 - xTemp2)/(currentFunct(xTemp2) - currentFunct(xTemp1))));
			
		xTemp1 = xTemp2;
		xTemp2 = xInter;
			
		*iterations = *iterations + 1;
	}
	while((fabs(currentFunct(xInter)) > tolerence) && (*iterations < maxIterations));
		
	return xInter;
	
}
/*	This function is the newton method, it returns the x value where a function crosses through zero
 *	First input is the initial guess for the newton method
 *	The second input is the tolerence for how close you need to be to beable to stop the method
 *	The third input is the pointer to an integer where the number of iterations the function goes through is stored
 *	The last input is an int representing the maximum number of iterations the method is allowed to go through to prevent infinite loops
 */
fpt newtonMethod (fpt x1, fpt tolerence, int* iterations, int maxIterations)
{
	fpt xNext = x1;
	
	*iterations = 0;
	
	do
	{
		*iterations = *iterations + 1;
		
		//find the next zero
		xNext = xNext - (currentFunct(xNext)/currentFunct_derivative(xNext));
		
	}
	while((fabs(currentFunct(xNext)) > tolerence) && (*iterations < maxIterations));
	
	return xNext;
}

int main()
{
	
	//set up variables
	fpt higherEnd = 2.1005;
	fpt lowerEnd = -1.1134;
	fpt tolerence = 0.000001;
	
	int iterations = 0;
	
	//test the bisection method
	fpt zero = bisectionMethod(lowerEnd, higherEnd, tolerence, &iterations, 100);
	printf("\tRoot from bisectionMethod is:\t%f,\t# of interations are: %d\n", zero, iterations);
	
	//test the regula-falsi method
	zero = regulaFalsiMethod(lowerEnd, higherEnd, tolerence, &iterations, 100);
	printf("\tRoot from regulaFalsiMethod is:\t%f,\t# of interations are: %d\n", zero, iterations);
	
	//test the secant method
	zero = secantMethod(lowerEnd, higherEnd, tolerence, &iterations, 100);
	printf("\tRoot from secantMethod is:\t%f,\t# of interations are: %d\n", zero, iterations);
	
	//test the newton method
	zero = newtonMethod(lowerEnd, tolerence, &iterations, 100);
	printf("\tRoot from newtonMethod is:\t%f,\t# of interations are: %d\n", zero, iterations);
	
	return 1;
}
 