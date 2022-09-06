% Function stacks different types of QCs in M_all



function M = func_stackingQCs(M_all, QC_names, QC_type, E)
	if strcmp(QC_names, 'all')
		QC_names = fieldnames(M_all);
	end

	switch QC_type
		case 'eQC'
			% initialize output
			M = [];

			for i = 1:length(QC_names)
				% check if QC_names{i} is in the M_all
				if isfield(M_all, QC_names{i})
					M = cat(3, M, M_all.(QC_names{i})(E));
				else
					warning([QC_names{i}, ' is not in M_all']);
				end
			end

		case 'ieQC'
			M.tilde = [];
			M.hat = [];
            M.L = [];

			for i = 1:length(QC_names)
				% check if QC_names{i} is in the M_all
				if isfield(M_all, QC_names{i})
					Mi = M_all.(QC_names{i})(E);

					M.tilde = cat(3, M.tilde, Mi.tilde);
					M.hat = cat(3, M.hat, Mi.hat);
                    
                    if isfield(Mi, 'L')
                        M.L = cat(3, M.L, Mi.L);
                    else
                        M.L = cat(3, M.L, []);
                    end
				else
					warning([QC_names{i}, ' is not in M_all']);
				end
			end
        
        otherwise
            error(['[QC_type] should be one of {EC,iEC}!']);
    end
end