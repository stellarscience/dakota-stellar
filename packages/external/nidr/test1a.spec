KEYWORD strategy {kwstart;0x2a;kwend;0x2b}
	[graphics {kwstart;0x33;kwend;0x34}]
	[ tabular_graphics_data {kwstart;0x23;kwend;0x45}
	  [tabular_graphics_file STRING {kwstart;0x123;kwend;0x456}] ]
	( multi_level {kwstart;0x24;kwend;0x25}
		( uncoupled {kwstart;9;kwend;10}
		  [ adaptive_hybrid {kwstart;3;kwend;4}
			progress_threshold REAL {kwstart;1;kwend;2} ]
                  [ num_solutions_transferred INTEGER {kwstart;5;kwend;6}]
		  method_list STRINGLIST {kwstart;7;kwend;8})
		|
		( coupled {kwstart;0x22;kwend;0x23}
		  global_method_pointer STRING{kwstart;0xb;kwend;0xc}
		  local_method_pointer STRING{kwstart;0xd;kwend;0xe}
		  [local_search_probability REAL{kwstart;0x20;kwend;0x21}] ) )
	|
	( surrogate_based_opt{kwstart;0x28;kwend;0x29}
		opt_method_pointer STRING {kwstart;0x26;kwend;0x27}
		[max_iterations INTEGER {kwstart;0x28;kwend;0x29}]
		)

KEYWORD method {kwstart;0x30;kwend;0x31}
	[id_method STRING {kwstart;0x2c;kwend;0x2d}]
	[model_pointer STRING {kwstart;0x2e;kwend;0x2f}]
;
