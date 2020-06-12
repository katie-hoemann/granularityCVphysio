function results = tablecompile_start(overallKey,newDataKey,newData,newDataCol,leftKeys,rightKeys,removeVar)
%joins iterative data to overall key variable, then deletes
%iterative key variable
results = horzcat(newDataKey,newData);
results = array2table(results,'VariableNames',{rightKeys,newDataCol});
results = outerjoin(overallKey,results,'LeftKeys',leftKeys,'RightKeys',rightKeys);
results = removevars(results,removeVar);
end