// Include the Arduino Stepper Library
#include <Stepper.h>

// Number of steps per output rotation
const int stepsPerRevolution = 100; // step distance will depend on your choice of motor. 

// Create Instance of Stepper library
Stepper myStepper(stepsPerRevolution, 8, 9, 10, 11);


void setup() 
// This loop is for a single stretch and relaxation cycle. Please feel free to adjust as needed for your own experiments.
// Potential options that cna be easily implemented include:
// Static strech 
{
  // set the speed at 60 rpm:
  myStepper.setSpeed(60);
  // initialize the serial port:
  Serial.begin(9600);

   Serial.println("clockwise"); // Turn the motor Clockwise (in the stretch position) 
   myStepper.step(stepsPerRevolution);
   delay(1000); // Hold delay for a few moments 

   Serial.println("counterclockwise"); // Turn the motor counterClockwise (Return to the unstretched position) 
   myStepper.step(-stepsPerRevolution);
   delay(1000); // Hold delay for a few moments 

   // set the speed at 60 rpm:
  myStepper.setSpeed(60);
  // initialize the serial port:
  Serial.begin(9600);

   Serial.println("clockwise");
   myStepper.step(stepsPerRevolution);
   delay(1000);
   Serial.println("counterclockwise");
   myStepper.step(-stepsPerRevolution);
   delay(1000);

}

void loop() 
