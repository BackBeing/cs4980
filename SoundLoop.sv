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

    parameter n = 16'b0000000000000001;
	
    // 3-bit register holding the current state of the FSM
    logic [2:0] state;
    
    // 3-bit combinational logic result
    logic [2:0] nextState;
	
	//holding the value of i and j
	logic [15:0] i;
	logic [3:0] j;
	
	//temp double holds
	logic [63:0] tempReal;
	logic [63:0] tempImag;
	logic [63:0] temp2Real;
	logic [63:0] temp2Imag;
    
    always_ff@(posedge Clock) begin
        if (Areset == 1'b0) begin
            state <= s00;
        end else begin
            state <= nextState;
        end
    end
    
    // Output logic
    // assign Out = (state == s11) ? 1'b1 : 1'b0;
    
    // Next state logic
    always_comb begin
        case(state)
        s000 : if (In == 1'b0) nextState = s000;
            else{
				i = 16'b1111111111111111;
				j = 4'b1111;
				nextState = s010;
			}
        s001 : nextState = s111;
        s010 : {
				i = 16'b1111111111111111;
				nextState = s100;
			}
        s011 : {j = j - 1;
				if (j == 0) nextState = s001;
				else nextState = s010;
				}
		s100 : {if (i & n) nextState = s101;
				else nextState = s110;
				}
        s101 : {i = i - 1;
				if (i == 0) nextState = s011;
				else nextState = s010;
				}
        s110 : {
	//	complex<double> temp = f[i];
	//		complex<double> Temp = W[(i * a) % (n * a)] * f[i + n];
	//		f[i] = temp + Temp;
	//		f[i + n] = temp - Temp;
				{
				tempReal = KgrpsReal [ 16'b1111111111111111 -i ];
				tempImag = KgrpsImag [ 16'b1111111111111111 -i ];
				temp2Real =  WReal [ ((16'b1111111111111111 -i )* 16'b0111111111111111) % (16'b0111111111111111* n ) ] * KgrpsReal[(16'b1111111111111111 -i ) + n];
				temp2Imag = WImag [ ((16'b1111111111111111 -i )* 16'b0111111111111111) % (16'b0111111111111111* n ) ] * KgrpsImag[(16'b1111111111111111 -i ) + n];
				KgrpsReal [ 16'b1111111111111111 -i ] = tempReal + temp2Real;
				KgrpsImag [ 16'b1111111111111111 -i ] = tempImag + temp2Imag;
				KgrpsReal[(16'b1111111111111111 -i ) + n] = tempReal - temp2Real;
				KgrpsImag[(16'b1111111111111111 -i ) + n] = tempImag - temp2Imag;
				
				}//all the actual computer steps
				//not assigned sign +/-
				
				i = i - 1;
				if (i == 0) nextState = s011;
				else nextState = s010;
				}
        s111 : if (In == 1'b0) nextState = s000;
				else nextState = s111;
        endcase
    end
endmodule
