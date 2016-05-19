function setSurfaceOverpotentialExpression(obj)
    switch obj.surfaceOverpotentialExpression
        case 'SternPotential'
            expression = '(phi_s-phi)*UT';
        case 'ZetaPotential'
            expression = '(phi_s-projectReactionPlaneToSurface(phi))*UT';
    end
    obj.surfaceVariables.set('V', expression, 'Metal - reaction plane potential difference, with dimension');
end