classdef    ValidatorFunctionCreator < handle     
%   Knows constraints and creates a validator function       
%
%   Keywords( ValidatorFunctionCreator )
%
%   See also: ValidatorFunctionCreator_test, InputPreprocessor_test

%   TODO:   compare the excecution time of ASSERT and VALIDATEATTRIBUTES with 
%           alternative contructs
%   TODO:   include an assert that validator is not {''}, i.e coloumn #5 
%
%   TODO:   Include new attribute values of the function, validateattributes.
%           e.g. 'decreasing','increasing','nondecreasing','nonincreasing' of R2016a
%
%   HISTORY
%   2011-03-21, poi: Handle truncated parameter names and values with **validatestring**
%   2011-07-10, poi: Added and removed 'IsoWeekNumber' to 'TimeUnitX'
%   2012-12-28, poi: MatlabVarName handles cell array of strings 
%   2012-12-28, poi: %#ok<DEFNU> as needed to make mlint button green
%   2013-02-12, poi: Use this function not only with InputputPreprocessor 
%   2013-07-30, poi: Added new type of constraint: USERSTRINGLISTS. 'member'
%   2013-07-30, poi: some refactoring. Now caValidator should be empty before the array 
%               of function_handles, fhv, is combined to one function_handle. 

%   HISTORY
%   2014-11-26, poi: fixed a silly copy&paste bug in function, member


%   Obsolete
%   lib.m       is replaced by ValidatorFunctionCreator  

%#ok<*AGROW>    allow variable to grow inside loop
%#ok<*NBRAK>    allow brackets, e.g. [1:12]

    properties( Constant = true, GetAccess = private )       

%       FUNCTION NAME LIST
        FunctionNameList = { 'CellNum', 'CellStr', 'Colormap', 'Datevec'    ...
                        ,    'FileExist', 'FileNotExist', 'FolderExist'     ...
                        ,    'FolderNotExist', 'MatlabVarName', 'HgHandle'  ...
                        ,    'addOptDummy'                                  }; 
                            
%       STRINGLISTS
%       --  These list of strings are used to validate string values
%       --  The names of the lists of strings are listed in the <StringLists> and then 
%           defined. Thus, a new name must appear twice. 
        StringLists = { 'TimeUnit', 'TimeUnitX' , 'Weekday' , 'TFI', 'PublishFormat'  ...
                    ,   'RGB'     , 'Permission','Encoding' , 'HgUnit'                };   

%       USERSTRINGLISTS
%       --  The list of strings, which is the second argument, is used to validate 
%           string values
%       -- 
        UserTwoArguments = { 'member' };
        
%       NUMERICVECTORS
%       --  These list of strings are used to validate numeric values
%       --  The names of the lists of strings are listed in the <NumericVectors> and  
%           then defined. Thus, a new name must appear twice. 
        NumericVectors  = { 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Weekday' };
        
%       Cell array of numerics, caNUMERIC
%       The names of these cell arrays of numerics are composed by a prefix, "ca", 
%       and a base name. The base name shall be a member of NumericVectors
        CellArrayNumeric    = { 'caMonth'   , 'caDay'   , 'caHour', 'caMinute'  ...
                            ,   'caSecond'  , 'caWeekday'                       }; 
    end
    
    properties( Constant = true, GetAccess = public )
%               
        ValidateAttributesClasses   =                       ...
            { 'numeric' , 'single'  , 'double'              ...
            , 'int8'    , 'int16'   , 'int32'   , 'int64'	... 
            , 'uint8'	, 'uint16'	, 'uint32'	, 'uint64'	...
            , 'logical'	, 'char'	, 'struct'	, 'cell'    ...
            , 'function_handle'     , 'class_name'	        };
        
        ValidateAttributesArguments =                               ...
            { '2d'      , 'integer'     , 'nonsparse'   , 'real'    ...
            , 'binary'                                              ...     R2009b
            , 'column'  , 'nonempty'    , 'nonzero'     , 'row'     ...   
            , 'even'    , 'nonnan'      , 'odd'         , 'scalar'  ...
            , 'finite'  , 'nonnegative' , 'positive'    , 'vector'  }; 
        
        ValidateAttributesTwoArguments =                                        ...
            { '<'   , '<='  , '>'   , '>='  , 'size', 'numel', 'ncols', 'nrows' };
        
%       The code rely on the *order* of the strings in the lists!
        TimeUnit    = { 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second' };
        
        TimeUnitX   =                                                           ...
                    { 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second'        ...
                    , 'IsoWeekday',                  'Holiday', 'Workingday'    ...
                    , 'UserTime', 'Analogue', 'Mode'                            };
                            
        IsoWeekday  = { 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun' };
        
        TFI         = { 'True', 'False', 'Ignore' };
        
        PublishFormat   = { 'html', 'ppt', 'xml', 'latex' }; %, 'doc', 'pdf' };
        
        RGB         = { 'Red', 'Green', 'Blue' }; 
        Permission  = { 'append', 'read', 'replace', 'write', 'r','w','a','r+','w+','a+'};
        
        Encoding    = { 'ibm850', 'latin1', 'latin9', 'US-ASCII', 'UTF-8' };
        
        HgUnit      = { 'pixels', 'normalized', 'inches', 'centimeters' ...
                    ,   'points', 'characters'                          };

        Month   = [1:1:12];  
        Day     = [1:1:31];
        Hour    = [0:1:23];
        Minute  = [0:1:59];
        Second  = [0:1:59];
        Weekday = [1:1: 7];     
        
    end  
    
    methods ( Access = public )     
        function    fh  = createValidatorFunction( this, caClass, caValidator ) 
%           Creates a function handle that either throws an exception or returns nothing.
%           This function handle is tailored for the method, parse, of the class,
%           inputParser. 
% 
%           VALIDATEATTRIBUTES is its first choice. Next is matching against a number 
%           of lists defined in this class, e.g. lists listed in this.StringLists. 
%          
            fhv = cell(0);
%           Firstly, convert those validator names, which refers to functions, to 
%           anonymous functions.

%           2013-07-29, poi: Added "isa( x, 'char' ) &&" to allow for constraints with
%           two arguments where the second is of any class
            ism = cellfun(  @(x)                                        ...
                            isa( x, 'char' )                            ...
                        &&  any( strcmp( x, this.FunctionNameList ) )   ...
                        ,   caValidator                                 ...
                        ,   'uni', true                                 );
                    
            for ca = caValidator( ism )
                fhv{1,end+1}    = @(x) feval( ca{1}, x );       
            end 
            caValidator = caValidator( not( ism ) );            % Remove  

%           Secondly, handle function handles to make it possible to use ISMEMBER       
            ism = cellfun( @(x) isa( x, 'function_handle' ), caValidator );
            for ca = caValidator( ism )
                fhv{1,end+1}    = @(x) assert( ca{1}(x)             ...
                    ,   'ValidatorFunctionCreator:FunctionHandle'   ...
                    ,   'Problem with function handle'              );
            end
            caValidator = caValidator( not( ism ) );    % Remove function handles
            
%           VALIDATEATTRIBUTES takes care of caClass and those elements of 
%           caValidator, which are members of ValidateAttributesArguments or
%           ValidateAttributesTwoArguments.

            is1 = cellfun( @(x)                                                 ...
                            isa( x, 'char' )                                    ...
                        &&  any( strcmp( x, this.ValidateAttributesArguments ) )...
                        ,   caValidator                                         ...
                        ,   'uni', true                                         );
                    
            is2 = cellfun( @(x)                                                     ...
                            isa( x, 'char' )                                        ...
                        &&  any( strcmp( x, this.ValidateAttributesTwoArguments ) ) ...
                        ,   caValidator                                             ...
                        ,   'uni', true                                             );
                    
            ism = is1 | is2 | [ false, is2(1:end-1) ];
                    
            if any( ism )
                fhv{1,end+1} = @(x) validateattributes( x, caClass, caValidator(ism));
            else
                fhv{1,end+1} = @(x) validateattributes( x, caClass, {} );
            end
            caValidator = caValidator( not( ism ) );
            
%           StringLists                  
%           ism = ismember( caValidator, this.StringLists );
            ism = cellfun( @(x)                                     ...
                            isa( x, 'char' )                        ...
                        &&  any( strcmp( x, this.StringLists ) )    ...
                        ,   caValidator                             ...
                        ,   'uni', true                             );
            
%           if any( ism )
            for ca = caValidator( ism )
                fhv{1,end+1} = @(x) assert( any( strcmp( x, this.(ca{1}) ) )    ...
                        , 'ValidatorFunctionCreator:UnknownName'                ...
                        , 'The name, "%s", is not a "%s"',  x, ca{1}            );
            end 
%           end
            caValidator = caValidator( not( ism ) );
           
%           NumericVectors            
%           ism = ismember( caValidator, this.NumericVectors );
%           ism = cellfun( @(x) any( strcmp( x, this.NumericVectors ) ), caValidator ); 
            ism = cellfun( @(x)                         ...
                isa( x, 'char' )                        ...
            &&  any( strcmp( x, this.NumericVectors ) ) ...
            ,   caValidator                             ...
            ,   'uni', true                             );

%           if any( ism )
            for ca = caValidator( ism )
                fhv{1,end+1} = @(x) assert( all( ismembc( x, this.(ca{1}) ) )    ...
                        , 'ValidatorFunctionCreator:UnknownNumber'               ...
                        , 'The number, "%i", is not a "%s"',  x, ca{1}           );
            end 
%           end
            caValidator = caValidator( not( ism ) );
        
%           Cell array of numerics, caNUMERIC
%           ism = ismember( caValidator, this.CellArrayNumeric );
%           ism = cellfun( @(x) any( strcmp( x, this.CellArrayNumeric ) ), caValidator );
            ism = cellfun( @(x)                             ...
                isa( x, 'char' )                            ...
            &&  any( strcmp( x, this.CellArrayNumeric ) )   ...
            ,   caValidator                                 ...
            ,   'uni', true                                 );
            
%           if any( ism )
            for ca = caValidator( ism )
                cacaStr         = regexp( ca{1}, '^ca(\w+)$', 'tokens' );
                NumVectorName   = cacaStr{1}{1};
                assert( ismember( NumVectorName, this.NumericVectors )      ...
                    ,   'ValidatorFunctionCreator:UnknownNumericVectorName' ...
                    ,   'The name, "%s", is not in the list, "%s"'          ...
                    ,   ca{1}, 'NumericVectors'                             )
                fhv{1,end+1} =                                                  ...
                    @(y) assert( all(                                           ...
                                cellfun( @(x)                                   ... 
                                        all( ismembc(x,this.(NumVectorName)) )  ...
                                        , y ) )                                 ...
                            , 'ValidatorFunctionCreator:UnknownNumberInCellArray' ...
                            , 'The caNumber, "%s", is not a "%s"'               ...
                            ,  '??? How to prints caNum ???', ca{1}             );
            end
            caValidator = caValidator( not( ism ) );
%           end
%
%           USERSTRINGLISTS
            ism = cellfun( @(x)                                         ...
                            isa( x, 'char' )                            ...
                        &&  any( strcmp( x, this.UserTwoArguments ) )   ...
                        ,   caValidator                                 ...
                        ,   'uni', true                                 );
           
            ixm = find( ism );
            for ix = ixm
                fhv{1,end+1} = @(x) feval( caValidator{ix:ix+1}, x );       
            end 
            caValidator([ixm,ixm+1]) = [];  % Remove  
            
%           Assert that all entries have been handled  
            assert( isempty( caValidator )                                          ...
            , 'ValidatorFunctionCreator:ParameterNotHandled'                        ...
            , ['Validator not handled: ', repmat('"%s", ', 1,length(caValidator) )] ...
            , caValidator{:}                                                        );
%{
            2013-07-30, poi: 
            assert( isempty( caValidator( not( ism ) ) )                            ...
            , 'ValidatorFunctionCreator:ParameterNotHandled'                        ...
            , ['Validator not handled: ', repmat('"%s", ', 1,numel(sum(not(ism))) )]...
            , caValidator{ not( ism ) }                                             );
%}                    
%        Combine all function handles of fhv        
%        fh  = @(x) cellfun( @(foo,v) feval(foo,v), fhv, repmat( {x}, 1,numel(fhv) ) );
%
%        The following construct, which is based on the function, CombineFunctionHandles, 
%        performs the same job as the previous cellfun-construct. The construct with 
%        CombineFunctionHandles is better. Creting the function_handle is somewhat faster 
%        with the cellfun-construct - a few tenth of a millisecond. However, evaluating
%        the function_handle is faster with the nested-function-construct. The difference
%        is significant when the function is called many times.

            fh  = @(x) CombineFunctionHandles( x );         
            function    CombineFunctionHandles( x )
                for fh_local = fhv
                    fh_local{:}( x )
                end
            end        
        end
        
        function    Keywords( this )
        %   Prints the keywords in the command window - simple reminder of the names
            
            clipboard( 'copy'                                           ,   ...
                sprintf(  [ 'persistent ipp\n'                              ...
                          , 'if nargin == 0,    return,     end\n'          ...
                          , 'if isempty( ipp )\n'                           ...
                          , '    ipp = InputPreprocessor( {\n'              ...
                          , '        1 ''Name''    ''''  {''char''}  {}\n'  ...
                          , '        } );\nend\n'                           ...
                          , 'ipv = ipp.parse( varargin{:} );\n' ]           ) )
               
            fprintf( '\nValidatorFunctionCreator, validateattributes\n' )
            for ca  = { 'ValidateAttributesClasses', 'ValidateAttributesArguments'  ...
                   ,    'ValidateAttributesTwoArguments'                            ...
                   ,    'FunctionNameList'                                          ...
                   ,    'StringLists', 'NumericVectors', 'CellArrayNumeric'         };
                
                fprintf( '%s', ca{:} )
                N  = numel( this.(ca{:}) );
                for ii = 1 : 5 : N
                    if ii <= N - 5
                        fprintf(  '\n    %-16s%-16s%-16s%-16s%-16s' ...
                                , this.(ca{:}){ii:ii+4}             )
                    else
                        fprintf(  ['\n    ', repmat('%-16s',1,N-ii+1)]  ...
                                , this.(ca{:}){ii:end}                  )
                    end
                end
                fprintf('\n')
            end          
        end
    end
end
%   -------------------------------------------------        
function    CellNum( caList )       %#ok<DEFNU>
    assert( all( cellfun( @(val) isnumeric(val), caList ) ) ...
        ,   'VFC:CellNum:NotCellNum'                        ...
        ,   'Not CellNum'                                   )
end
function    CellStr( caList )       %#ok<DEFNU>
    %   An extra ALL to handle matrices
    assert( all( all( cellfun( @(val) ischar(val), caList ) ) ) ...
        ,   'VFC:CellStr:NotCellStr'                            ...
        ,   'Not CellStr'                                       )
end
function    Colormap( vec )         %#ok<DEFNU>

    assert( size( vec, 2 ) == 3                             ...
        ,   'ValidatorFunctionCreator:Colormap:IllegalSize' ...
        ,   'Size of Colormap is "%s", but should be nx3'   ...
        ,   value2short( size( vec ) )                      )

    assert( max( vec(:) ) <= 1+eps                              ...
        ,   'ValidatorFunctionCreator:Colormap:IllegalValue'    ...
        ,   'Max value of Colormap is "%s", but must be <= 1'   )
    
    assert( min( vec(:) ) >= -eps                               ...
        ,   'ValidatorFunctionCreator:Colormap:IllegalValue'    ...
        ,   'Min value of Colormap is "%s", but must be >= 0'   )
    
end
%{
function    Datevec( vec )        
% TODO: vectorize, accept a column of Datavec
    assert( isa( vec, 'double' ), 'poi:isdatevec:IllegalType'       , ...
            'Inarg, "%s", should be a double', value2short( vec )   )
        
    assert( size(vec,1)==1 && ( numel(vec)==3 || numel(vec)==6 )                    ...
        ,   'poi:isdatevec:IllegalSize'                                             ...
        ,   'Size of inarg, "%s", should be 1x3 or 1x6', value2short( size(vec) )   )
    
    assert( isequal( vec, round( vec ) ), 'poi:isdatevec:NotFlint'  ...
        ,   'Inarg, "%s" should be flint', value2short( vec )       )
 
    assert( ismembc( vec(2), [ 1 : 12 ] ), 'poi:isdatevec:NotDatevec'   ...
        ,   'Second element of Inarg, "%s", should be [1:12]'           ...
        ,   value2short( vec )                                          )
    
    assert( ismembc( vec(3), [ 1 : 31 ] ), 'poi:isdatevec:NotDatevec'   ...
        ,   'Third element of Inarg, "%s", should be [1:31]'            ...
        ,   value2short( vec )                                          )
    
    if  numel(vec)==6

        assert( all( ismembc( vec(4), [ 0 : 23 ] ) ), 'poi:isdatevec:NotDatevec'    ...
            ,   'Two last elements of Inarg, "%s", should be [0:60]'                ...
            ,   value2short( vec )                                                  )
        
        assert( all( ismembc( vec(5:6), (0:60) ) ), 'poi:isdatevec:NotDatevec'      ...
            ,   'Two last elements of Inarg, "%s", should be [0:60]'                ...
            ,   value2short( vec )                                                  )
    end
end
%}
function    Datevec( vec )          %#ok<DEFNU>
%   Vectorized: accepts a column of Datevec

    assert( size(vec,2)==3 || size(vec,2)==6            ...
        ,   'vfc:isdatevec:IllegalSize'                 ...
        ,   'Size of inarg, "%s", should be nx3 or nx6' ...               
        ,   value2short( size(vec) )   )
    
    assert( isequal( vec, round( vec ) ), 'vfc:isdatevec:NotFlint'  ...
        ,   'Inarg, "%s" should be flint', value2short( vec )       )
 
    assert( all( ismembc( unique(vec(:,2)), [ 1 : 12 ] ) )      ...
        ,   'vfc:isdatevec:NotDatevec'                          ...
        ,   'Second element of Inarg, "%s", should be [1:12]'   ...
        ,   value2short( vec )                                  )
    
    assert( all( ismembc( unique(vec(:,3)), [ 1 : 31 ] ) )      ...
        ,   'vfc:isdatevec:NotDatevec'                          ...
        ,   'Third element of Inarg, "%s", should be [1:31]'    ...
        ,   value2short( vec )                                  )
    
    if  size(vec,2)==6

        assert( all( ismembc( unique(vec(:,4)), (0:23) ) )              ...
            ,   'vfc:isdatevec:NotDatevec'                              ...
            ,   'The forth element of Inarg, "%s", should be [0:23]'    ...
            ,   value2short( vec )                                      )
        
        assert( all( ismembc( unique(vec(:,5:6)), (0:60) ) )            ...
            ,   'vfc:isdatevec:NotDatevec'                              ...
            ,   'Two last elements of Inarg, "%s", should be [0:60]'    ...
            ,   value2short( vec )                                      )
    end
end
function    MatlabVarName( str )    %#ok<DEFNU>
    if iscellstr( str )
        isn = cellfun( @(s) isvarname(s), str, 'uni', true );
        assert( all( isn )                                  ...
            ,   'vfc:MatlabVarName:IllegalName'             ...
            ,   '"%s" is not a legal Matlab variable name'  ...
            ,   str( find( not( isn ), 1, 'first' ) )       )
    else
        assert( isvarname( str )                            ...
            ,   'vfc:MatlabVarName:IllegalName'             ...
            ,   '"%s" is not a legal Matlab variable name'  ...
            ,   str                                         )
    end
end
function    FileExist( str )        %#ok<DEFNU>
    assert( exist( str, 'file' ) == 2           ...
        ,   'vfc:FileExist2:CannotFindFile'     ...
        ,   'Cannot find file, "%s"'            ...
        ,   str                                 )
% FIXME: m-file without extension returns true, but fopen fails     
end
function    FileNotExist( str )     %#ok<DEFNU>
    assert( not( exist( str, 'file' ) == 2 )    ...
        ,   'vfc:FileExist2:FileFound'          ...
        ,   'File, "%s", exists'                ...
        ,   str                                 )
end
function    FolderExist( str )      %#ok<DEFNU>
    % 2014-07-16, poi: changed 'file' to 'dir'
    assert( exist( str, 'dir' ) == 7           ...
        ,   'vfc:FileExist7:CannotFindFolder'   ...
        ,   'Cannot find folder, "%s"'          ...
        ,   str                                 )
end
function    FolderNotExist( str )   %#ok<DEFNU>
    % 2014-07-16, poi: changed 'file' to 'dir' 
    assert( not( exist( str, 'dir' ) == 7 )    ...
        ,   'vfc:FileExist7:CannotFindFolder'   ...
        ,   'The folder, "%s" already exists'   ...
        ,   str                                 )
end
function    addOptDummy( ~ )        %#ok<DEFNU>
    assert( true, 'vfc:addOptDummy:Dummy', 'Dummy' )
end
function    HgHandle( num )         %#ok<DEFNU>
    assert( ishandle( num )                     ...
        ,   'vfc:HgHandle:WrongType'            ...
        ,   '"%s" isn''t handle graphic handle' ...
        ,   value2short( num )                  ) 
end
function    member( list, val )     %#ok<DEFNU>
%   assert( not( isempty( ismember( val, list ) ) )     ...     2014-11-26, poi: 
    assert( any( ismember( val, list ) )                ...
        ,   'vfc:member:NotFound'                       ...
        ,   '"%s" is not a member of the given list'    ...
        ,   value2short( val )                          )
end