`timescale 1ns / 1ps


module SoundLoop #(parameter WIDTH=64)(
    input [WIDTH - 1 :0] KgrpsReal [15:0],
	input [WIDTH - 1 :0] KgrpsImag [15:0],
	input [WIDTH - 1 :0] WReal [14:0],
	input [WIDTH - 1 :0] WImag [14:0],
	input In,
    input Clock,
    input Areset,
    output logic [16:0] i,
    output logic [4:0] j,
    output logic [WIDTH - 1 :0] KgrpsNewReal [15:0],
    output logic [WIDTH - 1 :0] KgrpsNewImag [15:0]
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
	
	/*
	reg [63:0] mem_kr [15:0];
    initial begin
    $readmemb("KgrpsRealDat", mem_2d);
    end
    reg [63:0] mem_ki [15:0];
    initial begin
    $readmemb("KgrpsImagDat", mem_2d2);
    end
    reg [63:0] mem_wr [14:0];
    initial begin
    $readmemb("WRealDat", mem_2d2);
    end
    reg [63:0] mem_wi [14:0];
    initial begin
    $readmemb("WImagDat", mem_2d2);
    end    
    */
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

	logic [16:0] nexti;

	logic [4:0] nextj;
	
	//temp double holds
	logic [WIDTH - 1:0] tempReal;
	logic [WIDTH - 1:0] tempImag;
	logic [WIDTH - 1:0] temp2Real;
	logic [WIDTH - 1:0] temp2Imag;
	
	logic [WIDTH * 2 :0] longdoubletp;
	logic signal;
    
    always_ff@(posedge Clock) begin
        if (Areset == 1'b0) begin
            state <= s000;
            i <= pp;
            j <= 4'b1111;
        end else begin
            state <= nextState;
            i <= nexti;
            j <= nextj;
        end
    end
    
    // Output logic
    // assign Out = (state == s11) ? 1'b1 : 1'b0;
    
    // Next state logic
    always_comb begin
        case(state)
        s000 : begin
                 nexti = pp;
                 nextj = 4'b1111; 
                 if (In == 1'b0) nextState = s000;
                 else nextState = s010;
               end
        s001 : begin
                nexti = i;
                nextj = j;
                nextState = s111;
               end
        s010 : begin
                nexti = pp;
                nextj = j;
                nextState = s100;
                end
        s011 : begin
                nexti = i;
                nextj = j - 1;
				if (nextj == 0) nextState = s001;
				else nextState = s010;
			    end
		s100 : begin 
		          nexti = i;
                  nextj = j;
                  qq = pp - i;
		          if (qq & n) nextState = s101;
				  else nextState = s110;
			     end
        s101 : begin
                nexti = nexti - 1;
                nextj = j;
				if (nexti == 0) nextState = s011;
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
				
				nexti = i - 1;
				nextj = j;
				if (nexti == 0) nextState = s011;
				else nextState = s010;
				end
        s111 : begin
                nexti = i;
                nextj = j;
                if (In == 1'b0) nextState = s000;
				else nextState = s111;
				end
        endcase
    end
endmodule
