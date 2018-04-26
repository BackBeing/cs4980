`timescale 1ns / 1ps


module SoundLoop(
    input [63:0] KgrpsReal [15:0],
	input [63:0] KgrpsImag [15:0],
	input [63:0] WReal [14:0],
	input [63:0] WImag [14:0],
	input In,
    output [63:0] KgrpsNewReal [15:0],
	output [63:0] KgrpsNewImag [15:0],
    input Clock,
    input Areset
    );
    
	//for(int j = 0; j < log2(N); j++) {
	//	for(int i = 0; i < N; i++) {
	//	  if(!(i & n)) {
	//		complex<double> temp = f[i];
	//		complex<double> Temp = W[(i * a) % (n * a)] * f[i + n];
	//		f[i] = temp + Temp;
	//		f[i + n] = temp - Temp;
	//	  }
	//	}
	//}
	
    // state encodings
    parameter s000 = 3'b000; //didn't enter j loop
	parameter s001 = 3'b001; //enter j loop and j = log2(N)
    parameter s010 = 3'b010; //enter j loop and j != log2(N)
	parameter s011 = 3'b011; //enter i loop and i = N
	parameter s100 = 3'b100; //enter i loop and i != N
    parameter s101 = 3'b101; //enter i loop but (i & n) = True
    parameter s110 = 3'b110; //enter i loop and !(i & n) 
	parameter s111 = 3'b111; //j loop done

    reg signed [15:0] n  = 16'b0000000000000001;
	reg signed [16:0] pp = 17'b01111111111111111;
	reg signed [15:0] kk = 16'b0111111111111111;
	
	logic [63:0] KgrpsNewReal [15:0];
    logic [63:0] KgrpsNewImag [15:0];
	
	
	logic [16:0] qq;
	logic [32:0] nn;
	logic [15:0] mm;
	logic [15:0] tp;
	logic [16:0] qres;
    // 3-bit register holding the current state of the FSM
    logic [2:0] state;
    
    // 3-bit combinational logic result
    logic [2:0] nextState;
	
	//holding the value of i and j
	logic [16:0] i;
	logic [4:0] j;
	
	//temp double holds
	logic [63:0] tempReal;
	logic [63:0] tempImag;
	logic [63:0] temp2Real;
	logic [63:0] temp2Imag;
	
	logic [126:0] longdoubletp;
	logic signal;
    
    always_ff@(posedge Clock) begin
        if (Areset == 1'b0) begin
            state <= s000;
        end else begin
            state <= nextState;
        end
    end
    
    // Output logic
    // assign Out = (state == s11) ? 1'b1 : 1'b0;
    
    // Next state logic
    always_comb begin
        case(state)
        s000 : if (In == 1'b0) begin
                    nextState = s000;
                    i = pp;
                    j = 4'b1111;
                    end
            else nextState = s010;
        s001 : nextState = s111;
        s010 : begin
                i = pp;
                nextState = s100;
                end
        s011 : begin
                j = j - 1;
				if (j == 0) nextState = s001;
				else nextState = s010;
			    end
		s100 : if (i & n) nextState = s101;
				else nextState = s110;
        s101 : begin
                i = i - 1;
				if (i == 0) nextState = s011;
				else nextState = s010;
				end
        s110 : begin
	//	complex<double> temp = f[i];
	//		complex<double> Temp = W[(i * a) % (n * a)] * f[i + n];
	//		f[i] = temp + Temp;
	//		f[i + n] = temp - Temp;
				//{
				qq = pp - i;
                nn = {qq[16],qq[15:0] * kk[14:0]};
                tp = (nn % kk);
                qres = qq + 1;
                
				tempReal = KgrpsReal [ qq ];
				tempImag = KgrpsImag [ qq ];
				
				//double 1 11 52
				longdoubletp = WReal [ tp[14:0] ][62:0] * KgrpsReal[ qres[15:0] ][62:0];
				signal = WReal[ tp[14:0] ][63] * KgrpsReal[ qres[15:0] ][63];
				temp2Real = {signal, longdoubletp[114:52] };
				
				longdoubletp = WImag [ tp[14:0] ][62:0] * KgrpsImag[ qres[15:0] ][62:0];
				signal = WImag [ tp[14:0] ][63] * KgrpsImag[ qres[15:0] ][63];
				temp2Imag = {signal, longdoubletp[114:52] };
				
				KgrpsNewReal [ qq ] = tempReal + temp2Real;
				KgrpsNewImag [ qq ] = tempImag + temp2Imag;
				KgrpsNewReal [ qres[15:0] ] = tempReal - temp2Real;
				KgrpsNewImag [ qres[15:0] ] = tempImag - temp2Imag;
				//}
				//all the actual computer steps
				//not assigned sign +/-
				
				i = i - 1;
				if (i == 0) nextState = s011;
				else nextState = s010;
				end
        s111 : if (In == 1'b0) nextState = s000;
				else nextState = s111;
        endcase
    end
endmodule
