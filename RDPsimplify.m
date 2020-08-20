function [ps,ix]=RDPsimplify(p, tol)

% Ramer-Douglas-Peucker algorithm for curve semplification

% [ps,ix] = dpsimplify(p,tol)
%
% dpsimplify uses the recursive Douglas-Peucker line simplification 
% algorithm to reduce the number of vertices in a piecewise linear curve 
% according to a specified tolerance. The algorithm is also know as
% Iterative Endpoint Fit. It works also for polylines and polygons
% in higher dimensions.
%
% In case of nans (missing vertex coordinates) dpsimplify assumes that 
% nans separate polylines. As such, dpsimplify treats each line
% separately.
%
% For additional information on the algorithm follow this link
% http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
%
% Input arguments
%
%     p     polyline n*d matrix with n vertices in d 
%           dimensions.
%     tol   tolerance (maximal euclidean distance allowed 
%           between the new line and a vertex)
%
% Output arguments
%
%     ps    simplified line
%     ix    linear index of the vertices retained in p (ps = p(ix))
%
% Examples
%
% 1. Simplify line 
%
%     tol    = 1;
%     x      = 1:0.1:8*pi;
%     y      = sin(x) + randn(size(x))*0.1;
%     p      = [x' y'];
%     ps     = dpsimplify(p,tol);
%
%     plot(p(:,1),p(:,2),'k')
%     hold on
%     plot(ps(:,1),ps(:,2),'r','LineWidth',2);
%     legend('original polyline','simplified')
%
% 2. Reduce polyline so that only knickpoints remain by 
%    choosing a very low tolerance
%
%     p = [(1:10)' [1 2 3 2 4 6 7 8 5 2]'];
%     p2 = dpsimplify(p,eps);
%     plot(p(:,1),p(:,2),'k+--')
%     hold on
%     plot(p2(:,1),p2(:,2),'ro','MarkerSize',10);
%     legend('original line','knickpoints')
%
% 3. Simplify a 3d-curve
% 
%     x = sin(1:0.01:20)'; 
%     y = cos(1:0.01:20)'; 
%     z = x.*y.*(1:0.01:20)';
%     ps = dpsimplify([x y z],0.1);
%     plot3(x,y,z);
%     hold on
%     plot3(ps(:,1),ps(:,2),ps(:,3),'k*-');
%
%
% 
% Author : Wolfgang Schwanghart, 13. July, 2010.
% Then modified by Roberto Rizzo @ HWU Edinburgh / University of Aberdeen,
% Date: April 2020

if nargin == 0
    help RPDsimplify
    return
end

error(nargchk(2, 2, nargin));

% error checking 
if ~isscalar(tol) || tol<0
    error('tollerance must be a positive scalar')
end

%nr of dimensions
nrvertices = size(p,1);
dims = size(p,1);

%anonymus function for starting-point and end-point comparison using a
%relative tollerance test
compare = @(a,b) abs (a-b)/max(abs(a), abs(b)) <= eps;

%what happens when we have NaNs?
%NaNs divide polylines.
Inan = any(isnan(p),2);
%any NaN at all?
Inanp = any(Inan);

%if there is only one vertex
if nrvertices == 1 || isempty(p)
    ps = p;
    ix = 1;
    
%if there are two vertices
elseif nrvertices == 2 && ~Inanp
    %when the line has no vertices (except end- and starting-point of the line)
    %check if the distance between both is less that the tolerance.
    %If so, return the centre.
    if dims == 2
        d = hypot(p(1,1)-p(2,1),p(1,2)-p(2,2));
    else
        d = sqrt(sum((p(1,:)-p(2,:)).^2));
    end
    
    if d <= tol
        ps = sum(p,1)/2;
        ix = 1;
    else
        ps = p;
        ix = [1,2];
    end
    
elseif Inanp
    
    %case: there are NaNs in the p array
    %--> find start and end indices of contiguous non-NaNs data
    Inan = ~Inan;
    sIX = strfind(Inan', [0,1])' + 1;
    eIX = strfind(Inan', [1,0])';
    
    if Inan(end) == true
        eIX = [eIX;nrvertices];
    end
    
    if Inan(1)
        sIX = [1;sIX];
    end
    
    % calculate lenght on non-NaNs components
    lIX = eIX - sIX+1;
    %put each component in a single cell
    c = mat2cell(p(Inan,:), lIX,dims);
    
    %now call RDPsimplyfy again inside cellfun.
    if nargout == 2
        [ps,ix] = cellfun(@(x) RDPsimplify(x,tol),c,'UniformOutput', false);
        ix = cellfun(@(x,six) x+six-1,ix, num2cell(sIX),'UniformOutput', false);
    else
        ps = cellfun(@(x) RDPsimplify(x,tol),c,'UniformOutput',false);
    end
    
    %write data for the cell array back to a matrix
    ps = cellfun(@(x) [x;nan(1,dims)], ps, 'UniformOutput', false);
    ps = cell2mat(ps);
    ps(end,:) = [];
    
    %write ix into a matrix
    if nargout == 2;
        ix = cell2mat(ix);
    end
    
    
else
    
    %if there are no NaNs then start the recursive algorithm
    
    ixe = size(p,1);
    ixs =  1;
    
    %logical vector for the vertices to be retained
    I = true(ixe,1);
    
    %call recursive function
    p = simplifyrec(p,tol,ixs,ixe);
    ps = p(I,:);
    
    %if desired return the index of retained vertices
    if nargout == 2;
        ix = find(I);
    end
end

% ================== RECURSIVE FUNCTION ===============================

function p = simplifyrec(p,tol,ixs,ixe)
    
    % check if starting-point and end-point are the same better comparison
    % needed which included tolerance eps
    
    c1 = num2cell(p(ixs,:));
    c2 = num2cell(p(ixe,:));
    
    %same start- and end-point with tolerance
    sameSE = all(cell2mat(cellfun(compare,c1(:),c2(:),'UniformOutput',false)));
    
    if sameSE
        %calculate the shortest distance of all vertices between ixs and
        %ixe to ixs only
        if dims == 2
            d = hypot(p(ixs,1)-p(ixs+1:ixe-1,1),p(ixs,2)-p(ixs+1:ixe-1,2));
        else
            d    = sqrt(sum(bsxfun(@minus,p(ixs,:),p(ixs+1:ixe-1,:)).^2,2));
        end
    else
        % calculate the shortest distance of all points to the line from
        % ixs to ixe; subtract starting point from other locations.
        pt = bsxfun(@minus, p(ixs+1:ixe,:), p(ixs,:));
        
        % end-point
        a = pt(end,:)';
        
        beta = (a' * pt')./(a'*a);
        b = pt-bsxfun(@times,beta,a)';
        
        if dims == 2
            % if line in 2D use the numerical more robust hypot function
            d = hypot(b(:,1),b(:,2));
        else
            d = sqrt(sum(b.^2,2));
        end
    end
    
    % identify maximum distance and get the linear index of its location
    [dmax, ixc] = max(d);
    ixc = ixs +ixc;
    
    % if the max distance is smaller the the tollerance remore vertices
    % between ixs and ixe
    if dmax <= tol
        if ixs ~= ixe-1
            I(ixs+1:ixe-1) = false;
        end
    % if not, call symplifyrec for the segment between ixs and ixc (ixc and ixe)
    else
        p = simplifyrec(p,tol,ixs,ixc);
        p = simplifyrec(p,tol,ixc,ixe);
    end
end
end

    