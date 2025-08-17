/*
  GA25-370 + Quadrature Encoder + L298N (Minimal ID Logger)
  CSV: t,duty,Vs,omega,theta   (omega = RPM, theta = revolutions)

  Commands @115200:
    SEQ        -> run 20 s step test (0–5s:0%, 5–15s:50%, 15–20s:0%)
    STOP       -> stop sequence (hold last duty)
    Dxx        -> set duty, e.g., D60 or D-40 (exits sequence)
    CSV n      -> set sample period in ms (10..1000), e.g., CSV 10
    CPR n      -> set counts-per-rev for your ISR/wiring, e.g., CPR 120
*/

//////////////// Pins ////////////////
const int ENA = 5;   // PWM
const int IN1 = 8;
const int IN2 = 7;

const int ENC_A = 2; // must be interrupt-capable pin
const int ENC_B = 4;

///////////// Config /////////////////
float SUPPLY_VOLTAGE   = 12.0;   // printed in CSV
float COUNTS_PER_REV   = 120.0;  // 1x decoding (A rising only)
int   ENCODER_POLARITY = -1;     // flip sign if needed
unsigned long SAMPLE_INTERVAL_MS = 10;

struct Step { float t_s; int duty; };
Step seq[] = { {0.0,0}, {5.0,80}, {15.0,0} };
const float SEQ_T_END = 20.0;

////////// State //////////
volatile long encCount = 0;
unsigned long t0_ms=0, lastSampleMs=0;
long lastCount=0;
int  currentDutyCmd = 0;
bool runSequence = false;

////////// Motor //////////
void setMotorPercent(int dutyPercent){
  int d = constrain(abs(dutyPercent),0,100);
  int pwm = map(d,0,100,0,255);
  if (dutyPercent >= 0){ digitalWrite(IN1,HIGH); digitalWrite(IN2,LOW); }
  else                { digitalWrite(IN1,LOW);  digitalWrite(IN2,HIGH); }
  analogWrite(ENA,pwm);
}

////////// Encoder ISR (A rising) //////////
void onEncA(){
  bool b = digitalRead(ENC_B);
  int step = (b==LOW)? +1 : -1;
  encCount += ENCODER_POLARITY * step;
}

////////// Serial commands //////////
void handleSerial(){
  if(!Serial.available()) return;
  String cmd = Serial.readStringUntil('\n'); cmd.trim();
  if(!cmd.length()) return;

  if(cmd.startsWith("D")||cmd.startsWith("d")){
    int val = cmd.substring(1).toInt();
    currentDutyCmd = constrain(val,-100,100);
    runSequence = false;
    setMotorPercent(currentDutyCmd);
    Serial.print("# duty="); Serial.println(currentDutyCmd);

  }else if(cmd.equalsIgnoreCase("SEQ")){
    runSequence = true; t0_ms = millis();
    Serial.println("# sequence start");

  }else if(cmd.equalsIgnoreCase("STOP")){
    runSequence = false; setMotorPercent(currentDutyCmd);
    Serial.println("# sequence stop");

  }else if(cmd.startsWith("CSV")){
    long ms = cmd.substring(3).toInt();
    if(ms>=10 && ms<=1000){ SAMPLE_INTERVAL_MS = (unsigned long)ms; Serial.print("# CSV="); Serial.println((long)SAMPLE_INTERVAL_MS); }

  }else if(cmd.startsWith("CPR")){
    long cpr = cmd.substring(3).toInt();
    if(cpr>0 && cpr<50000){ COUNTS_PER_REV = (float)cpr; Serial.print("# CPR="); Serial.println((long)COUNTS_PER_REV); }
  }
}

////////// Setup/Loop //////////
void setup(){
  pinMode(ENA,OUTPUT); pinMode(IN1,OUTPUT); pinMode(IN2,OUTPUT);
  pinMode(ENC_A,INPUT_PULLUP); pinMode(ENC_B,INPUT_PULLUP);
  Serial.begin(115200);
  attachInterrupt(digitalPinToInterrupt(ENC_A), onEncA, RISING);
  setMotorPercent(0);
  t0_ms = millis(); lastSampleMs = t0_ms;

  // one-time header (MATLAB logger skips it if it appears mid-stream)
  Serial.println("t,duty,Vs,omega,theta");
}

void loop(){
  handleSerial();

  unsigned long now = millis();
  float t = (now - t0_ms)/1000.0f;

  // sequence profile
  if(runSequence){
    int duty = seq[0].duty;
    if(t >= seq[1].t_s) duty = seq[1].duty;
    if(t >= seq[2].t_s) duty = seq[2].duty;
    setMotorPercent(duty);
    if(t >= SEQ_T_END){ runSequence=false; setMotorPercent(0); Serial.println("# sequence done"); }
  }

  // periodic CSV
  if(now - lastSampleMs >= SAMPLE_INTERVAL_MS){
    unsigned long dt_ms = now - lastSampleMs; lastSampleMs = now;

    noInterrupts(); long cnt = encCount; interrupts();
    long dCounts = cnt - lastCount; lastCount = cnt;

    double dRev   = (double)dCounts / COUNTS_PER_REV;
    double theta  = (double)cnt / COUNTS_PER_REV;  // revolutions
    double dt     = dt_ms / 1000.0;
    double rpm    = dRev * 60.0 / dt;              // RPM

    int dutyOut = runSequence ? (t<5?0:(t<15?50:0)) : currentDutyCmd;

    float tPrint = (now - t0_ms)/1000.0f;
    Serial.print(tPrint,3); Serial.print(',');
    Serial.print(dutyOut);  Serial.print(',');
    Serial.print(SUPPLY_VOLTAGE,2); Serial.print(',');
    Serial.print(rpm,6);    Serial.print(',');
    Serial.println(theta,6);
  }
}