function Jp = jacobian_hexagoniser(J,Mesh)
% J is the Jacobian on the tetrahedral mesh, Mesh is the return value of
% minecraft(Mesh,d)
    Jp = zeros(size(J,1),size(Mesh.Hex,1));
    for i = 1:size(Jp,2)
        %Jp(:,i) = sum(J(:,Mesh.cells{i}.p),2);
        Jp(:,i) = sum(J(:,Mesh.cells{i}),2);
    end