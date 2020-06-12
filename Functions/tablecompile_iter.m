function results = tablecompile_iter(oldData,newDataKey,newData,newDataCol,leftKeys,rightKeys,removeVar)
%joins iterative data to overall key variable, then deletes
%iterative key variable
results2 = horzcat(newDataKey,newData);
results2 = array2table(results2,'VariableNames',{rightKeys,newDataCol});
results = outerjoin(oldData,results2,'LeftKeys',leftKeys,'RightKeys',rightKeys);
results = removevars(results,removeVar);
end